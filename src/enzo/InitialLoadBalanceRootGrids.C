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
void icol(int *x, int n, int m, FILE *log_fptr);
void fpcol(Eflt64 *x, int n, int m, FILE *fptr);
int DetermineNumberOfNodes(void);

int InitialLoadBalanceRootGrids(FILE *fptr, int TopGridRank,
				int TopGridDim, int &NumberOfRootGrids,
				int* &RootProcessors)
{

  if (NumberOfProcessors == 1 || LoadBalancing <= 1)
    return SUCCESS;

  /* Declarations */

  bool FinishedLevelZero;
  char line[MAX_LINE_LENGTH];
  int i, dim, GridID, Rank, ThisLevel, dummy, GridDims[MAX_DIMENSION];
  int ThisTask, *GridMap;
  int size, grid_num, RootGridID, GridPosition[MAX_DIMENSION];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  double *NumberOfSubgridCells;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ThisCellWidth;

  /* Determine number of nodes and thus cores/node */

  int NumberOfNodes = DetermineNumberOfNodes();

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

  // Compute (processor number)->(grid number) map
  Enzo_Dims_create(NumberOfRootGrids, TopGridRank, LayoutTemp);
  for (dim = 0; dim < TopGridRank; dim++)
    Layout[TopGridRank-1-dim] = LayoutTemp[dim];
  for (dim = TopGridRank; dim < MAX_DIMENSION; dim++)
    Layout[TopGridRank-1-dim] = 0;

  GridMap = new int[NumberOfRootGrids];
  RootProcessors = new int[NumberOfRootGrids];
  NumberOfSubgridCells = new double[NumberOfRootGrids];
  for (i = 0; i < NumberOfRootGrids; i++)
    NumberOfSubgridCells[i] = 0;

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

      ThisCellWidth = (RightEdge[0] - LeftEdge[0]) /
       (GridDims[0] - 2*DEFAULT_GHOST_ZONES);
      ThisLevel = nint(-logf(TopGridDim * ThisCellWidth) / M_LN2);

      // Record GridMap number if level-0
      if (ThisLevel == 0) {

       for (dim = 0; dim < Rank; dim++)
         GridPosition[dim] =
           int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
               (DomainRightEdge[dim] - DomainLeftEdge[dim]));
       grid_num = GridPosition[0] +
         Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
       GridMap[grid_num] = GridID-1;
       RootProcessors[GridID-1] = ThisTask;

      } // ENDIF level 0

      // From the grid map, determine in which root grid the subgrid
      // is located and add to subgrid cell count.
      if (ThisLevel == 1) {

       for (dim = 0, size = 1; dim < Rank; dim++) {
         GridPosition[dim] =
           int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
               (DomainRightEdge[dim] - DomainLeftEdge[dim]));
         size *= GridDims[dim];
       }
       grid_num = GridPosition[0] +
         Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
       RootGridID = GridMap[grid_num];
       NumberOfSubgridCells[RootGridID] += size;

      } // ENDIF level 1

    } // ENDIF GravityBoundaryType read

  } // ENDWHILE lines

  /* Now that we have the number of subgrid cells in each topgrid, we
     can distribute the topgrid tiles so each compute node have a
     similar total number of subgrid cells.  Place them cyclical and
     then move them around. */

  int Node, NewProc;
  double *CellsPerNode = new double[NumberOfNodes];
  int *GridsPerNode = new int[NumberOfNodes];
  double TotalSubgridCells = 0;

  for (i = 0; i < NumberOfNodes; i++) {
    CellsPerNode[i] = 0;
    GridsPerNode[i] = 0;
  }

  for (i = 0; i < NumberOfRootGrids; i++) {
    if (LoadBalancing == 2)
      Node = RootProcessors[i] / CoresPerNode;
    else if (LoadBalancing == 3)
      Node = RootProcessors[i] % NumberOfNodes;
    CellsPerNode[Node] += NumberOfSubgridCells[i];
    TotalSubgridCells += NumberOfSubgridCells[i];
    GridsPerNode[Node]++;
  }

  // Load balance the nodes
  bool Done = false;
  int MinNode, MaxNode;
  double MaxVal, MinVal;
  int MaxGridsPerNode = NumberOfRootGrids / NumberOfNodes;
  double DesiredSubgridCells = TotalSubgridCells / NumberOfNodes;
  double OptimalSize, DifferenceFromOptimal;

  while (!Done) {

    MaxNode = -1;
    MinNode = -1;
    MaxVal = 0;
    MinVal = huge_number;

    for (i = 0; i < NumberOfNodes; i++) {
      if (CellsPerNode[i] > MaxVal) {// || GridsPerNode[i] > MaxGridsPerNode) {
       MaxVal = CellsPerNode[i];
       MaxNode = i;
      }
      if (CellsPerNode[i] < MinVal) {// && GridsPerNode[i] < MaxGridsPerNode) {
       MinVal = CellsPerNode[i];
       MinNode = i;
      }
    }

//    if (debug) {
//      printf("LBRoot: Min/Max Node = %d(%lg)/%d(%lg)\n",
//            MinNode, MinVal, MaxNode, MaxVal);
//      printf("Cells:");
//      fpcol(CellsPerNode, NumberOfNodes, 16, stdout);
//      printf("Grids:");
//      icol(GridsPerNode, NumberOfNodes, 16, stdout);
//    }

    if (MaxVal > LOAD_BALANCE_RATIO*MinVal) {
      for (i = 0; i < NumberOfRootGrids; i++) {

       if (LoadBalancing == 2)
         Node = RootProcessors[i] / CoresPerNode;
       else if (LoadBalancing == 3)
         Node = RootProcessors[i] % NumberOfNodes;

       OptimalSize = (DesiredSubgridCells - MinVal) / 
	 (MaxGridsPerNode - GridsPerNode[MinNode]);
       DifferenceFromOptimal = fabs(OptimalSize - NumberOfSubgridCells[i]) /
	 OptimalSize;

       if (Node == MaxNode &&
           (NumberOfSubgridCells[i] < (MaxVal-MinVal)/2 ||
	    DifferenceFromOptimal < 0.5)) {

//	    (GridsPerNode[MinNode] < MaxGridsPerNode &&
//	     GridsPerNode[MaxNode] > MaxGridsPerNode))) {

         if (LoadBalancing == 2)  // block scheduling
           RootProcessors[i] = MinNode * CoresPerNode + GridsPerNode[MinNode];
         else if (LoadBalancing == 3)  // round-robin
           RootProcessors[i] = MinNode + NumberOfNodes * GridsPerNode[MinNode];

         /* Update node subgrid cell counts */

         CellsPerNode[MinNode] += NumberOfSubgridCells[i];
         CellsPerNode[MaxNode] -= NumberOfSubgridCells[i];

	 GridsPerNode[MinNode]++;
	 GridsPerNode[MaxNode]--;

         break;

       } // ENDIF move
      } // ENDFOR root grids

      // If we didn't find an appropriate transfer then quit.
      Done = (i == NumberOfRootGrids);

    } else
      Done = true;

  } // ENDWHILE !Done

  rewind(fptr);
  delete [] GridMap;
  delete [] NumberOfSubgridCells;
  delete [] CellsPerNode;
  delete [] GridsPerNode;

  return SUCCESS;
}
