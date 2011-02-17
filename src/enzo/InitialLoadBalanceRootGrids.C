/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE ROOT GRIDS
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include <string.h>
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
#include "communication.h"
#include "CommunicationUtilities.h"

#define LOAD_BALANCE_RATIO 1.05

int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int DetermineNumberOfNodes(void);
int LoadBalanceSimulatedAnnealing(int NumberOfGrids, int NumberOfNodes, 
				  int* &ProcessorNumbers, double* NumberOfCells, 
				  double* NumberOfSubcells);
int LoadBalanceHilbertCurveRootGrids(FLOAT *GridCenters[], int *CellCount,
				     int NumberOfGrids, int* &RootProcessors);

int InitialLoadBalanceRootGrids(FILE *fptr, int TopGridRank,
				int TopGridDim, int &NumberOfRootGrids,
				int* &RootProcessors)
{

  /* Determine number of nodes and thus cores/node */

  int NumberOfNodes = DetermineNumberOfNodes();

  if (NumberOfProcessors == 1 || LoadBalancing == 0)
    return SUCCESS;

  /* If this is a zoom-in calculation, we shouldn't load balance the
     root grids by nodes but the finest static grid. */

//  if (StaticRefineRegionLevel[0] != INT_UNDEFINED)
//    return SUCCESS;

  /* Declarations */

  bool FinishedLevelZero;
  char line[MAX_LINE_LENGTH];
  int i, dim, GridID, Rank, ThisLevel, dummy, GridDims[MAX_DIMENSION];
  int ThisTask, *GridMap;
  int *Workload;
  int size, grid_num, RootGridID, GridPosition[MAX_DIMENSION];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  double *NumberOfCells, *NumberOfSubgridCells;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT *GridCenters[MAX_DIMENSION];
  float ThisCellWidth;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Layout[dim] = 0;
    GridPosition[dim] = 0;
  }

  // Root processor does all of the work, and broadcasts the new
  // processor numbers

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  // Count number of level-0 grids (assume that all level-0 grids at
  // the beginning)
  NumberOfRootGrids = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    if (sscanf(line, "Grid = %"ISYM, &dummy) > 0)
      NumberOfRootGrids++;

    // Reached the last root grid
    if (strstr(line, "NextGridNextLevel") != NULL)
      break;

  } // ENDWHILE lines

  } // ENDIF root processor

  // If we're doing normal load balancing, synchronize and exit.
  if (LoadBalancing == 1) {
    rewind(fptr);
#ifdef USE_MPI
    MPI_Bcast(&NumberOfRootGrids, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif /* USE_MPI */
    return SUCCESS;
  } // ENDIF LoadBalancing == 1

  if ((LoadBalancing == 2 || LoadBalancing == 3) &&
      MyProcessorNumber == ROOT_PROCESSOR) {

  // Compute (processor number)->(grid number) map
  Enzo_Dims_create(NumberOfRootGrids, TopGridRank, LayoutTemp);
  for (dim = 0; dim < TopGridRank; dim++)
    Layout[TopGridRank-1-dim] = LayoutTemp[dim];
  for (dim = TopGridRank; dim < MAX_DIMENSION; dim++)
    Layout[TopGridRank-1-dim] = 0;

  RootProcessors = new int[NumberOfRootGrids];
  GridMap = new int[NumberOfRootGrids];
  NumberOfCells = new double[NumberOfRootGrids];
  NumberOfSubgridCells = new double[NumberOfRootGrids];
  for (i = 0; i < NumberOfRootGrids; i++) {
    NumberOfCells[i] = 0;
    NumberOfSubgridCells[i] = 0;
  }

  // Start reading hierarchy
  rewind(fptr);
  FinishedLevelZero = false;
  Rank = 1;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    sscanf(line, "Grid = %"ISYM, &GridID);
    sscanf(line, "Task = %"ISYM, &ThisTask);
    sscanf(line, "GridRank = %"ISYM, &Rank);
    switch (Rank) {
    case 1:
      sscanf(line, "GridDimension = %"ISYM, GridDims);
      sscanf(line, "GridLeftEdge = %"PSYM, LeftEdge);
      sscanf(line, "GridRightEdge = %"PSYM, RightEdge);
      break;
    case 2:
      sscanf(line, "GridDimension = %"ISYM" %"ISYM, GridDims, GridDims+1);
      sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM, LeftEdge, LeftEdge+1);
      sscanf(line, "GridRightEdge = %"PSYM" %"PSYM, RightEdge, RightEdge+1);
      break;
    case 3:
      sscanf(line, "GridDimension = %"ISYM" %"ISYM" %"ISYM,
            GridDims, GridDims+1, GridDims+2);
      sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM" %"PSYM,
            LeftEdge, LeftEdge+1, LeftEdge+2);
      sscanf(line, "GridRightEdge = %"PSYM" %"PSYM" %"PSYM,
            RightEdge, RightEdge+1, RightEdge+2);
      break;
    default:
      ENZO_FAIL("Bad grid rank!");
    } // ENDSWITCH Rank

    // Compute level and other derived quantities after reading the
    // last line
    if (sscanf(line, "GravityBoundaryType = %"ISYM, &dummy) > 0) {

      for (dim = 0, size = 1; dim < Rank; dim++) {
	GridPosition[dim] =
	  int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
	      (DomainRightEdge[dim] - DomainLeftEdge[dim]));
	size *= GridDims[dim];
      }
      grid_num = GridPosition[0] +
	Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);

      ThisCellWidth = (RightEdge[0] - LeftEdge[0]) /
       (GridDims[0] - 2*DEFAULT_GHOST_ZONES);
      ThisLevel = nint(-logf(TopGridDim * ThisCellWidth) / M_LN2);

      // Record GridMap number if level-0
      if (ThisLevel == 0) {

       GridMap[grid_num] = GridID-1;
       RootProcessors[GridID-1] = ThisTask;
       NumberOfCells[GridID-1] += size;

      } // ENDIF level 0

      // From the grid map, determine in which root grid the subgrid
      // is located and add to subgrid cell count.
      else if (ThisLevel == 1) {

       RootGridID = GridMap[grid_num];
       NumberOfSubgridCells[RootGridID] += size;

      } // ENDIF level 1

    } // ENDIF GravityBoundaryType read

  } // ENDWHILE lines

  /* Now that we have the number of subgrid cells in each topgrid, we
     can distribute the topgrid tiles so each compute node have a
     similar total number of subgrid cells.  Use simulated annealing
     to load balance. */

  LoadBalanceSimulatedAnnealing(NumberOfRootGrids, NumberOfNodes, 
				RootProcessors, NumberOfCells,
				NumberOfSubgridCells);

  delete [] GridMap;
  delete [] NumberOfCells;
  delete [] NumberOfSubgridCells;

  } // ENDIF ROOT_PROCESSOR and LoadBalancing == 2 or 3

  if (MyProcessorNumber == ROOT_PROCESSOR && LoadBalancing == 4) {
    RootProcessors = new int[NumberOfRootGrids];
    Workload = new int[NumberOfRootGrids];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      GridCenters[dim] = new FLOAT[NumberOfRootGrids];

    for (i = 0; i < NumberOfRootGrids; i++)
      Workload[i] = 0;

    // Start reading hierarchy
    rewind(fptr);
    Rank = 1;
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

      sscanf(line, "Grid = %"ISYM, &GridID);
      sscanf(line, "Task = %"ISYM, &ThisTask);
      sscanf(line, "GridRank = %"ISYM, &Rank);
      switch (Rank) {
      case 1:
	sscanf(line, "GridDimension = %"ISYM, GridDims);
	sscanf(line, "GridLeftEdge = %"PSYM, LeftEdge);
	sscanf(line, "GridRightEdge = %"PSYM, RightEdge);
	break;
      case 2:
	sscanf(line, "GridDimension = %"ISYM" %"ISYM, GridDims, GridDims+1);
	sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM, LeftEdge, LeftEdge+1);
	sscanf(line, "GridRightEdge = %"PSYM" %"PSYM, RightEdge, RightEdge+1);
	break;
      case 3:
	sscanf(line, "GridDimension = %"ISYM" %"ISYM" %"ISYM,
	       GridDims, GridDims+1, GridDims+2);
	sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM" %"PSYM,
	       LeftEdge, LeftEdge+1, LeftEdge+2);
	sscanf(line, "GridRightEdge = %"PSYM" %"PSYM" %"PSYM,
	       RightEdge, RightEdge+1, RightEdge+2);
	break;
      default:
	ENZO_FAIL("Bad grid rank!");
      } // ENDSWITCH Rank

      // Compute level and other derived quantities after reading the
      // last line
      if (sscanf(line, "GravityBoundaryType = %"ISYM, &dummy) > 0) {

	for (dim = 0, size = 1; dim < Rank; dim++) {
	  GridCenters[dim][GridID-1] =
	    (LeftEdge[dim] - DomainLeftEdge[dim]) /
	    (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	  size *= GridDims[dim];
	}

	Workload[GridID-1] = size;

	if (GridID == NumberOfRootGrids)
	  break;

      } // ENDIF GravityBoundaryType read

    } // ENDWHILE lines

    LoadBalanceHilbertCurveRootGrids(GridCenters, Workload, 
				     NumberOfRootGrids, RootProcessors);

    delete [] Workload;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] GridCenters[dim];

  } // ENDIF ROOT_PROCESSOR and LoadBalancing == 4

#ifdef USE_MPI
  MPI_Bcast(&NumberOfRootGrids, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
  if (NumberOfRootGrids > 1) {
    if (MyProcessorNumber != ROOT_PROCESSOR)
      RootProcessors = new int[NumberOfRootGrids];
    MPI_Bcast(RootProcessors, NumberOfRootGrids, IntDataType, ROOT_PROCESSOR, 
	      MPI_COMM_WORLD);
  } else {
    delete [] RootProcessors;
    RootProcessors = NULL;
  }
#endif

  rewind(fptr);

  return SUCCESS;
}
