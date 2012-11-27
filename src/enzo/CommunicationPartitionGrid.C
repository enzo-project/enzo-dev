/***********************************************************************
/
/  COMMUNICATION ROUTINE: PARTITION GRID
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       March, 2006
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
 
// Function prototypes
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int LoadBalanceHilbertCurve(grid *GridPointers[], int NumberOfGrids, 
			    int* &NewProcessorNumber);

#define USE_OLD_CPU_DISTRIBUTION

/* This option code ensures that in nested grid sims, the children grids are not split between two grids at level-1.
   It is off by default. */
#define CONTAINED_WITHIN_PARENT_OFF

#ifdef CONTAINED_WITHIN_PARENT
int *AllStartIndex[MAX_STATIC_REGIONS][MAX_DIMENSION];
int AllDims[MAX_STATIC_REGIONS][MAX_DIMENSION];
FLOAT AllLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
bool FirstTimeCalled = true;
#endif

int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum)
{
 
  if (NumberOfProcessors*NumberOfRootGridTilesPerDimensionPerProcessor == 1)
    return SUCCESS;
 
  // Declarations
 
  int Rank, dim, i, j, k, ijk, Dims[MAX_DIMENSION], Layout[] = {0,0,0};
  int TempDims[MAX_DIMENSION] /* , TempStart[MAX_DIMENSION] */ ;
  int *GridDims[MAX_DIMENSION], *StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
 
  if (debug) printf("Enter CommunicationPartitionGrid.\n");

  /* Initialize storage for grid left edges */

  CommunicationBarrier();
 
  /* Attach RandomForcingFields as BaryonFields (for the duration
     of partitioning only). */
 
  if (RandomForcing == 1 && ParallelRootGridIO != 1)
    Grid->GridData->AppendForcingToBaryonFields(); //AK
 
  /* Compute side length of root grid. */
 
  Grid->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
 
  for (dim = 0; dim < Rank; dim++)
    Dims[dim] -= 2*NumberOfGhostZones;
 
  float Edge = POW(float(Dims[0]*Dims[1]*Dims[2])/float(NumberOfProcessors),
		   1/float(Rank));
 
 
  /* If using MPI, use their routine to calculate layout. */
 
#ifdef USE_MPI
 
  int LayoutTemp[] = {0,0,0};

/*
  MPI_Arg Nnodes = NumberOfProcessors;
  MPI_Arg Ndims = Rank;
  MPI_Arg LayoutDims[] = {0, 0, 0};
 
  if (MPI_Dims_create(Nnodes, Ndims, LayoutDims) != MPI_SUCCESS) {
    ENZO_FAIL("Error in MPI_Dims_create.");
  }
*/

  int Nnodes = NumberOfProcessors;
  int Ndims = Rank;
  int LayoutDims[] = {0, 0, 0};

  if (Enzo_Dims_create(Nnodes, Ndims, LayoutDims) != SUCCESS) {
    ENZO_FAIL("Error in Enzo_Dims_create.");
  }

  for (dim = 0; dim < Rank; dim++)
    LayoutTemp[dim] = LayoutDims[dim];
 
  /* Swap layout because we want smallest value to be at Layout[0]. */
 
  for (dim = 0; dim < Rank; dim++)
    Layout[dim] = LayoutTemp[Rank-1-dim] * NumberOfRootGridTilesPerDimensionPerProcessor;
 
  /* Force some distributions if the default is brain-dead. */

/*
  if (Rank == 3 && NumberOfProcessors == 8)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 2;

  if (Rank == 3 && NumberOfProcessors == 64)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 4;

  if (Rank == 3 && NumberOfProcessors == 125)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 5;
 
  if (Rank == 3 && NumberOfProcessors == 216)
    for (dim = 0; dim < Rank; dim++)
      Layout[dim] = 6;
*/

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", Layout[0], Layout[1], Layout[2]);
  }

#endif /* USE_MPI */
 
 
  /* Generate arrays of grid dimensions and start positions. */
 
  int NumberOfNewGrids = 1;
  int DisplacementCount;
  float ExactDims, ExactCount;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
 
    /* Compute number of new grids along this dimension. */
 
    if (Layout[dim] == 0)
      Layout[dim] = max(nint((float)(Dims[dim])/Edge), 1);
 
    GridDims[dim] = new int[Layout[dim]];
    StartIndex[dim] = new int[Layout[dim]];
 
    /* Compute dims and start indexes of the grids along this dim. */
 
    ExactDims = float(Dims[dim])/float(Layout[dim]);
    ExactCount = 0;
    DisplacementCount = 0;
 
    for (i = 0; i < Layout[dim]; i++) {
      ExactCount += ExactDims;
 
      /* Compute the integer number of cells along this dimension
	 (if dim == 0 then make sure it is even as well since the FFT
	 requires this). */
 
      if (dim == 0)
	GridDims[dim][i] = nint(ExactCount*0.5)*2 - DisplacementCount;
      else
	GridDims[dim][i] = nint(ExactCount) - DisplacementCount;
 
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += GridDims[dim][i];
    }
 
    NumberOfNewGrids *= Layout[dim];
 
  }

  /* If specified, ensure that the new grids are contained within the
     partitions of the parent grid */

#ifdef CONTAINED_WITHIN_PARENT

  int ParentGridNum, NumberOfCoarseSlabs, CoarseSlab, ThisSlab, FirstCoarseSlab;
  int ThisStartIndex, ThisEndIndex, ThisWidth, NumberOfSlabs, CoarseLimit;
  int NumberOfCoarseSlabsLeft, ThisLevel;
  int ParentDims[MAX_DIMENSION];
  int *CoarseEdges;
  float ExactDimsLeft, ThisExactDims;
  FLOAT ParentLeftEdge[MAX_DIMENSION], ParentRightEdge[MAX_DIMENSION];
  bool FoundIt;
  HierarchyEntry *TempGrid;

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  if (FirstTimeCalled) {
    for (i = 0; i < MAX_STATIC_REGIONS; i++) {
      for (j = 0; j < MAX_DIMENSION; j++) {
	AllStartIndex[i][j] = NULL;
	AllLeftEdge[i][j] = FLOAT_UNDEFINED;
	AllDims[i][j] = INT_UNDEFINED;
      }
    }
    FirstTimeCalled = false;
  }

  // Determine this level by going up the AMR tree
  ThisLevel = 0;
  TempGrid = Grid;
  while (TempGrid->ParentGrid != NULL) {
    TempGrid = TempGrid->ParentGrid;
    ThisLevel++;
  }

  if (Grid->ParentGrid != NULL) {

    // Find parent in static regions.  If level==1, then it has to be
    // the root grid.
    ParentGridNum = INT_UNDEFINED;
    
    if (ThisLevel == 1) {
      ParentGridNum = 0;
    } else {
    
      Grid->ParentGrid->GridData->ReturnGridInfo(&Rank, ParentDims, ParentLeftEdge, 
						 ParentRightEdge);
      for (i = 0; i < MAX_STATIC_REGIONS; i++) {
	// If we're at the end of the list, stop
	if (StaticRefineRegionLevel[i] == INT_UNDEFINED) break;
	FoundIt = true;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  FoundIt &= (ParentLeftEdge[dim] >= StaticRefineRegionLeftEdge[i][dim]) &&
	    (ParentRightEdge[dim] <= StaticRefineRegionRightEdge[i][dim]) &&
	    (ThisLevel-1 == StaticRefineRegionLevel[i]+1);
	if (FoundIt) {
	  ParentGridNum = i+1;
	  break;
	}
      } // ENDFOR static regions

    } // ENDELSE ThisLevel == 1

    if (ParentGridNum == INT_UNDEFINED) {
      ENZO_VFAIL("CommunicationPartitionGrid: grid %d (%d), Parent not found?\n",
	      gridnum, ThisLevel)
    }

    for (dim = 0; dim < MAX_DIMENSION; dim++) {

      // If there's only one slab, it's easy
      if (Layout[dim] == 1) {
	GridDims[dim][0] = Dims[dim];
	StartIndex[dim][0] = 0;
	continue;
      }

      CoarseEdges = new int[Layout[dim]];

      // Calculate partitioned parent edges in units of this grid's
      // cell width.  Assume that gridnum = 0 is the root partitioned
      // grid.
      for (i = 0; i < Layout[dim]; i++)
	CoarseEdges[i] = 2 * AllStartIndex[ParentGridNum][dim][i] -
	  nint((Left[dim] - AllLeftEdge[ParentGridNum][dim]) * 
	       AllDims[0][dim] * POW(RefineBy, ThisLevel));

      // Count the number of coarse slabs that overlap with this grid
      NumberOfCoarseSlabs = 0;
      for (i = 0; i < Layout[dim]-1; i++)
	if ((CoarseEdges[i] >= 0 && CoarseEdges[i] < Dims[dim]) ||
	    (CoarseEdges[i+1] >= 0 && CoarseEdges[i+1] < Dims[dim]))
	  NumberOfCoarseSlabs++;
      if (CoarseEdges[Layout[dim]-1] >= 0 && CoarseEdges[Layout[dim]-1] < Dims[dim])
	NumberOfCoarseSlabs++;

      // Determine the first coarse slab that overlaps
      CoarseSlab = 0;
      ThisSlab = 0;
      while (CoarseEdges[CoarseSlab] <= 0 && CoarseSlab+1 < Layout[dim])
	CoarseSlab++;
      FirstCoarseSlab = CoarseSlab;

      /* Loop over all contained coarse slabs and split them
	 accordingly.  There must exist some buffer between grids,
	 i.e. original child and parent grids cannot have the same
	 edges. */

      while (ThisSlab < Layout[dim]) {

	// Get left and right edge of coarse slab, making sure that
	// they don't go outside of the original grid.
	if (CoarseSlab < NumberOfCoarseSlabs+FirstCoarseSlab &&
	    CoarseSlab < Layout[dim])
	  ThisEndIndex = min(CoarseEdges[CoarseSlab], Dims[dim]);
	else
	  ThisEndIndex = Dims[dim];
	ThisStartIndex = max(CoarseEdges[CoarseSlab-1], 0);
	ThisWidth = ThisEndIndex - ThisStartIndex;

	// How many coarse slabs are left
	NumberOfCoarseSlabsLeft = NumberOfCoarseSlabs - (CoarseSlab - FirstCoarseSlab);
	
	// Optimal splitting size and number of blocks within this coarse slab
	ExactDimsLeft = float(Dims[dim] - ThisStartIndex) / (Layout[dim] - ThisSlab);
	NumberOfSlabs = max( nint( float(ThisWidth) / ExactDimsLeft ), 1);
	  
	// Ensure that we don't over-split the current coarse slab so
	// we have enough slabs left for the remaining coarse slabs.
	CoarseLimit = NumberOfCoarseSlabsLeft - (Layout[dim] - ThisSlab) + 1;
	if (CoarseLimit > 0 && ThisEndIndex != Dims[dim])
	  NumberOfSlabs = min(NumberOfSlabs, CoarseLimit);
	
	// Now we can split the coarse slab
	ThisExactDims = float(ThisWidth) / NumberOfSlabs;
	ExactCount = 0.0;
	DisplacementCount = 0;

	for (i = 0; i < NumberOfSlabs; i++) {
	  ExactCount += ThisExactDims;
	  GridDims[dim][ThisSlab+i] = nint(0.5*ExactCount)*2 - DisplacementCount;
	  StartIndex[dim][ThisSlab+i] = ThisStartIndex + DisplacementCount;
	  DisplacementCount += GridDims[dim][ThisSlab+i];
	} // ENDFOR split coarse slab
	
	CoarseSlab++;
	ThisSlab += NumberOfSlabs;

      } // ENDWHILE slabs

      delete [] CoarseEdges;

    } // ENDFOR dims

  } // ENDIF parent exists

  /* Store this level's new partitions for use later */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    AllStartIndex[gridnum][dim] = new int[Layout[dim]];
    AllLeftEdge[gridnum][dim] = Left[dim];
    AllDims[gridnum][dim] = Dims[dim];
    for (i = 0; i < Layout[dim]; i++)
      AllStartIndex[gridnum][dim][i] = StartIndex[dim][i];
  }

  } // ENDIF ROOT_PROCESSOR

  // Broadcast values to all processors
#ifdef USE_MPI
  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    MPI_Bcast(GridDims[dim], Layout[dim], DataTypeInt, ROOT_PROCESSOR, MPI_COMM_WORLD);
    MPI_Bcast(StartIndex[dim], Layout[dim], DataTypeInt, ROOT_PROCESSOR, MPI_COMM_WORLD);
  }
#endif  

#endif /* CONTAINED_WITHIN_PARENT */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
    printf("PartitionGrid (on all processors): Layout = %"ISYM" %"ISYM" %"ISYM"\n",
      Layout[0], Layout[1], Layout[2]);
    printf("NumberOfNewGrids = %"ISYM"\n",NumberOfNewGrids);
 
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
      printf("GridDims[%"ISYM"]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %"ISYM,GridDims[dim][i]);
      }
      printf("\n");
    }
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
      printf("StartIndex[%"ISYM"]: ",dim);
      for (i = 0; i < Layout[dim]; i++)
      {
        printf(" %"ISYM,StartIndex[dim][i]);
      }
      printf("\n");
    }
    }

/*
  if ((ProblemType == 30) && (ParallelRootGridIO == 1) && (ParallelParticleIO == 1))
  {
    printf("Unigrid: %"ISYM"\n", Unigrid);
    printf("Set Unigrid = 1\n");
    Unigrid = 1;
  }
*/
 
  /* Initialize the under subgrid field for particle movement. */
 
  if (!ParallelRootGridIO) {
    if (debug) printf("Call ZeroSUS on TopGrid\n");
    Grid->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  }
 
  /* Generate this many grids (on this processor). */
 
/*
  Unigrid = 0;
  if (debug) printf("Re-set Unigrid = 0\n");
*/
 
  if (debug) printf("Grid structure: %"ISYM"\n", (int) (sizeof(grid)));
  if (debug) printf("SubGrids structure: %"ISYM"\n", (int) ((Layout[0]*Layout[1]*Layout[2])*sizeof(grid)));
 
  grid *NewGrid, *OldGrid = Grid->GridData;
  grid **SubGrids = new grid*[Layout[0]*Layout[1]*Layout[2]];
  HierarchyEntry *ThisGrid;
 
  int gridcounter = 0;
 
  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {
 
	/* Allocate a new grid hierarchy entry and insert into linked list. */
 
	if (gridcounter == 0)
	  ThisGrid = Grid;
	else {
	  ThisGrid = new HierarchyEntry;
	  ThisGrid->NextGridThisLevel = Grid->NextGridThisLevel;
	  ThisGrid->NextGridNextLevel = NULL;
	  ThisGrid->ParentGrid        = Grid->ParentGrid;
	  Grid->NextGridThisLevel     = ThisGrid;
	}
 
	/* Allocate a new grid and prepare it. */
 
	NewGrid = new grid;
	ThisGrid->GridData = NewGrid;
	NewGrid->InheritProperties(OldGrid);
	NewGrid->SetGravityParameters(OldGrid->ReturnGravityBoundaryType());
 
	/* Compute grid region. */
 
//      printf("GC K J I: %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",gridcounter,k,j,i);
 
	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  ijk = (dim == 0) ? i : ((dim == 1) ? j : k);
	  TempDims[dim] = GridDims[dim][ijk];
	  LeftEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk])/FLOAT(Dims[dim]);
	  RightEdge[dim] = Left[dim] + (Right[dim] - Left[dim])*
	    FLOAT(StartIndex[dim][ijk]+TempDims[dim])/FLOAT(Dims[dim]);
	  if (dim < Rank)
	    TempDims[dim] += 2*NumberOfGhostZones;
 
//        printf("  LeftEdge[%"ISYM"] = %8.4"FSYM"  RightEdge[%"ISYM"] = %8.4"FSYM"\n",
//               dim, LeftEdge[dim], dim, RightEdge[dim]);
 
	}
 
	NewGrid->PrepareGrid(Rank, TempDims, LeftEdge, RightEdge, 0);
 
	/* Record this subgrid number in the oldgrid's undersubgrid field. */
 
//      printf("Call ZeroSUS on OldGrid with Value = %10.4e\n", float(gridcounter+1));
	if (!ParallelRootGridIO)
	  if (OldGrid->ZeroSolutionUnderSubgrid(NewGrid,
		   ZERO_UNDER_SUBGRID_FIELD, float(gridcounter+1)) == FAIL) {
	    ENZO_FAIL("Error in grid->ZeroSolutionUnderSubgrid.");
	  }
	SubGrids[gridcounter] = NewGrid;
 
	gridcounter++;
 
      }
 
 
  Unigrid = 0;
  if (debug) printf("Re-set Unigrid = 0\n");
 
  /* Move Particles (while still on same processor). */

  if (!ParallelRootGridIO)
    if (OldGrid->MoveSubgridParticlesFast(gridcounter, SubGrids, TRUE) == FAIL) {
      ENZO_FAIL("Error in grid->MoveSubgridParticlesFast.");
    }
 
  int *PartitionProcessorNumbers = NULL;
  if (LoadBalancing == 4)
    LoadBalanceHilbertCurve(SubGrids, gridcounter,
			    PartitionProcessorNumbers);

  delete [] SubGrids;

  /* Distribute new grids amoung processors (and copy out fields). */

  CommunicationBarrier();
 
  gridcounter = 0;
  ThisGrid = Grid;

  if (MyProcessorNumber == ROOT_PROCESSOR) 
    printf("Grid distribution\n");
 
  for (k = 0; k < Layout[2]; k++)
    for (j = 0; j < Layout[1]; j++)
      for (i = 0; i < Layout[0]; i++) {
 
	grid *NewGrid = ThisGrid->GridData;
 
	/* Broadcast the number of particles to the other processors
	   (OldGrid is assumed to be on the root processor). */
 
	if (NumberOfProcessors > 1) {

	  int IntTemp = NewGrid->ReturnNumberOfParticles();
 
//          printf("NewGrid->ReturnNumberOfParticles: %"ISYM"\n", IntTemp);
 
	  CommunicationBroadcastValue(&IntTemp, ROOT_PROCESSOR);

	  NewGrid->SetNumberOfParticles(IntTemp);

//          printf("NG particle number set to %"ISYM"\n", IntTemp);

	}
 
	/* Transfer from Old to New (which is still also on root processor) */
 
	FLOAT Zero[] = {0,0,0};
 
        if (ParallelRootGridIO == FALSE) {
	  if (MyProcessorNumber == ROOT_PROCESSOR)
	    NewGrid->AllocateGrids();
 
          if (NewGrid->CopyZonesFromGrid(OldGrid, Zero) == FAIL) {
            ENZO_FAIL("Error in grid->CopyZonesFromGrid.");
          }

        } // ENDIF no PartitionNestedGrids
 
 
	/* Set processor number of new grid.  Cyclic distribution. */
 
        int NewProc = gridcounter % NumberOfProcessors;
        int ProcMap = ABS(NewProc - NumberOfProcessors) % NumberOfProcessors;
 
        if(NewGrid->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge) == FAIL) {
          ENZO_FAIL("Error in grid->ReturnGridInfo.");
        }

	if (PartitionProcessorNumbers != NULL)
	  NewProc = PartitionProcessorNumbers[gridcounter];
 
	/* Move Grid from current processor to new Processor. */
 
#ifdef USE_OLD_CPU_DISTRIBUTION
	if (!ParallelRootGridIO)
	  NewGrid->CommunicationMoveGrid(NewProc);
	else
	  NewGrid->SetProcessorNumber(NewProc);
#endif

#ifdef USE_PERMUTED_CPU_DISTRIBUTION
	if (!ParallelRootGridIO)
	  NewGrid->CommunicationMoveGrid(ProcMap);
	else
	  NewGrid->SetProcessorNumber(ProcMap);
#endif

	// TA introduced to be able to restart from large unigrid files
	// but allow parallel rootgrid IO
#ifdef USE_NEW_CPU_DISTRIBUTION
	  NewGrid->CommunicationMoveGrid(NewProc);
#endif
	

	// some debug output
        if (MyProcessorNumber == ROOT_PROCESSOR && debug1) {
          printf("Grid = %"ISYM", K J I: [%"ISYM",%"ISYM",%"ISYM"] Proc = %"ISYM"\n", gridcounter, k, j, i, NewProc);
          for (dim = 0; dim < Rank; dim++) {
            printf("  %"ISYM" ::  LeftEdge[%"ISYM"] = %8.4"PSYM"  RightEdge[%"ISYM"] = %8.4"PSYM"\n",
                   NewProc, dim, LeftEdge[dim], dim, RightEdge[dim]);
          }
	}

        /* Detach ForcingFields from BaryonFields. */
 
        if (RandomForcing == 1 && ParallelRootGridIO != 1)
          NewGrid->RemoveForcingFromBaryonFields(); //AK
 
	gridcounter++;
	ThisGrid = ThisGrid->NextGridThisLevel;
 
      }

  CommunicationBarrier();
 
  /* Clean up. */

  if (RandomForcing == 1 && ParallelRootGridIO != 1)
    OldGrid->RemoveForcingFromBaryonFields();
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Delete OldGrid\n");
 
  delete OldGrid;
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("OldGrid deleted\n");

  delete [] PartitionProcessorNumbers;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] GridDims[dim];
    delete [] StartIndex[dim];
  }

  if (debug) printf("Exit CommunicationPartitionGrid.\n");

  CommunicationBarrier();
 
  return SUCCESS;
}
