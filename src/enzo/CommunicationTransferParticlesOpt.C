/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PARTICLES
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified:   May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
 
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
#include "SortCompareFunctions.h"
void my_exit(int status);
 
// function prototypes
 
Eint32 compare_grid(const void *a, const void *b);
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int CommunicationShareParticles(int *NumberToMove, particle_data* &SendList,
				int &NumberOfReceives,
				particle_data* &SharedList);
int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);

#define NO_DEBUG_CTP
#define KEEP_PARTICLES_LOCAL
 
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids,
				   int TopGridDims[])
{

  if (NumberOfGrids == 1)
    return SUCCESS;
 
  int i, j, jstart, jend, dim, grid, proc, DisplacementCount, ThisCount;
  float ExactDims, ExactCount;

  /* Assume that the grid was split by Enzo_Dims_create, and create a
     map from grid number to an index that is determined from (i,j,k)
     of the grid partitions. */

  int *GridMap = new int[NumberOfGrids];
  int *StartIndex[MAX_DIMENSION];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  int Rank, grid_num, CenterIndex, bin;
  int *pbin;
  int GridPosition[MAX_DIMENSION], Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Layout[dim] = 0;
    GridPosition[dim] = 0;
  }

  GridPointer[0]->ReturnGridInfo(&Rank, Dims, Left, Right); // need rank
  if (Enzo_Dims_create(NumberOfGrids, Rank, LayoutTemp) == FAIL) {
    ENZO_FAIL("Error in Enzo_Dims_create.\n");
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Layout[dim] = 0;
  for (dim = 0; dim < Rank; dim++)
    Layout[Rank-1-dim] = LayoutTemp[dim];

  /* For unequal splits of the topgrid, we need the start indices of
     the partitions.  It's easier to recalculate than to search for
     them. */

  for (dim = 0; dim < Rank; dim++) {

    StartIndex[dim] = new int[Layout[dim]+1];
    ExactDims = float(TopGridDims[dim]) / float(Layout[dim]);
    ExactCount = 0.0;
    DisplacementCount = 0;

    for (i = 0; i < Layout[dim]; i++) {
      ExactCount += ExactDims;
      if (dim == 0)
	ThisCount = nint(0.5*ExactCount)*2 - DisplacementCount;
      else
	ThisCount = nint(ExactCount) - DisplacementCount;
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += ThisCount;
    } // ENDFOR i

    StartIndex[dim][Layout[dim]] = TopGridDims[dim];

  } // ENDFOR dim

  for (grid = 0; grid < NumberOfGrids; grid++) {
    GridPointer[grid]->ReturnGridInfo(&Rank, Dims, Left, Right);
    for (dim = 0; dim < Rank; dim++) {

      if (Layout[dim] == 1) {
	GridPosition[dim] = 0;
      } else {

	CenterIndex = 
	  int(TopGridDims[dim] *
	      (0.5*(Right[dim]+Left[dim]) - DomainLeftEdge[dim]) /
	      (DomainRightEdge[dim] - DomainLeftEdge[dim]));
      
	GridPosition[dim] = 
	  search_lower_bound(StartIndex[dim], CenterIndex, 0, Layout[dim],
			     Layout[dim]);

      } // ENDELSE

    } // ENDFOR dim
//    if (debug)
//      printf("grid %d: GPos = %d %d %d\n", grid, GridPosition[0],
//	     GridPosition[1], GridPosition[2]);
    grid_num = GridPosition[0] + 
      Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
    GridMap[grid_num] = grid;
  } // ENDFOR grids
 
  int *NumberToMove = new int[NumberOfProcessors];
  particle_data *SendList = NULL;
  particle_data *SharedList = NULL;
 
  for (i = 0; i < NumberOfProcessors; i++)
    NumberToMove[i] = 0;

#ifdef DEBUG_CTP
  if (MyProcessorNumber == ROOT_PROCESSOR) 
  for (j = 0; j < NumberOfGrids; j++)
    fprintf(stderr, "Pa(%"ISYM") CTP grid[%"ISYM"] = %"ISYM"\n",
	    MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfParticles());
#endif
 
  /* Generate the list of particle moves. */

  int Zero = 0;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->
	CommunicationTransferParticles(GridPointer, NumberOfGrids, grid, 
				       TopGridDims, NumberToMove, Zero, 
				       Zero, SendList, Layout, StartIndex, 
				       GridMap, COPY_OUT) == FAIL) {
      ENZO_FAIL("Error in grid->CommunicationTransferParticles(OUT).\n");
    }

  int TotalNumberToMove = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];

  int NumberOfReceives = 0;
#ifdef KEEP_PARTICLES_LOCAL

  /* Sort by grid number.  We no longer share with other processors
     because the particles are stored locally until
     CommunicationCollectParticles(ALLGRIDS). */

  SharedList = SendList;
  NumberOfReceives = TotalNumberToMove;
  int particle_data_size = sizeof(particle_data);
  //qsort(SharedList, TotalNumberToMove, particle_data_size, compare_grid);
  std::sort(SharedList, SharedList+TotalNumberToMove, cmp_grid());

#else

  /* Share the Particle moves. */

  if (CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
				  SharedList) == FAIL) {
    ENZO_FAIL("Error in CommunicationShareParticles.\n");
  }

#endif

  /* Copy particles back to grids */

  jstart = 0;
  jend = 0;

  // Copy shared particles to grids, if any
  if (NumberOfReceives > 0) {
    for (j = 0; j < NumberOfGrids && jend < NumberOfReceives; j++) {
      while (SharedList[jend].grid <= j) {
	jend++;
	if (jend == NumberOfReceives) break;
      }
      if (GridPointer[j]->CommunicationTransferParticles
	  (GridPointer, NumberOfGrids, j, TopGridDims,
	   NumberToMove, jstart, jend, 
	   SharedList, Layout, StartIndex, GridMap, COPY_IN) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationTransferParticles(IN).\n");
      }
      jstart = jend;
    } // ENDFOR grids
  } // ENDIF NumberOfRecieves > 0

  /* Even if we have no receives, we still have to remove the sent
     particles. */

#ifndef KEEP_PARTICLES_LOCAL
  else {

    for (j = 0; j < NumberOfGrids; j++)
      if (GridPointer[j]->ReturnProcessorNumber() == MyProcessorNumber)
	if (GridPointer[j]->CleanUpMovedParticles() == FAIL) {
	  ENZO_FAIL("Error in grid->CleanUpMovedParticles.\n");
	}

  } // ENDELSE NumberOfReceives > 0
#else

  /* If the loop with grid::CTP stopped before we covered all of the
     grids, be sure to remove the particles on the remaining grids on
     all processors. If NumberOfReceives == 0, then no particles were
     moved. */

  if (NumberOfReceives > 0)
    for (j = SharedList[NumberOfReceives-1].grid; j < NumberOfGrids; j++)
      if (GridPointer[j]->CleanUpMovedParticles() == FAIL) {
	ENZO_FAIL("Error in grid->CleanUpMovedParticles.\n");
      }

#endif /* KEEP_PARTICLES_LOCAL */
 
  /* Set number of particles so everybody agrees. */

#ifdef UNUSED
  if (NumberOfProcessors > 1) {
    int *AllNumberOfParticles = new int[NumberOfGrids];
    for (j = 0; j < NumberOfGrids; j++)
      if (GridPointer[j]->ReturnProcessorNumber() == MyProcessorNumber)
	AllNumberOfParticles[j] = GridPointer[j]->ReturnNumberOfParticles();
      else
	AllNumberOfParticles[j] = 0;

    CommunicationAllSumValues(AllNumberOfParticles, NumberOfGrids);
    for (j = 0; j < NumberOfGrids; j++)
      GridPointer[j]->SetNumberOfParticles(AllNumberOfParticles[j]);

    delete [] AllNumberOfParticles;
  }
#endif /* UNUSED */

  /* Cleanup. */

  if (SendList != SharedList)
    delete [] SendList;
  delete [] SharedList;
  delete [] NumberToMove;
  delete [] GridMap;
  for (dim = 0; dim < Rank; dim++)
    delete [] StartIndex[dim];

  CommunicationSumValues(&TotalNumberToMove, 1);
  if (debug)

    printf("CommunicationTransferParticles: moved = %"ISYM"\n",
  	   TotalNumberToMove);
 
  return SUCCESS;
}

