/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE ROOT GRIDS BY SUBGRID CELLS
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

int CommunicationLoadBalanceRootGrids(LevelHierarchyEntry *LevelArray[], 
				      int TopGridRank, int CycleNumber)
{

  if (NumberOfProcessors == 1 || LoadBalancing <= 1)
    return SUCCESS;

  if (CycleNumber % LoadBalancingCycleSkip != 0)
    return SUCCESS;

  /* Declarations */

  LevelHierarchyEntry *Temp;

  char line[MAX_LINE_LENGTH];
  int i, dim, GridID, Rank, ThisLevel, dummy, ThisTask, GridDims[MAX_DIMENSION];
  int size, grid_num, RootGridID, NumberOfRootGrids;
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION], GridPosition[MAX_DIMENSION];
  double *NumberOfSubgridCells;
  int *RootProcessors, *GridMap;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  // Count number of level-0 grids
  NumberOfRootGrids = 0;
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
    NumberOfRootGrids++;

  // Allocate some memory
  GridMap = new int[NumberOfRootGrids];
  RootProcessors = new int[NumberOfRootGrids];
  NumberOfSubgridCells = new double[NumberOfRootGrids];

  for (i = 0; i < NumberOfRootGrids; i++)
    NumberOfSubgridCells[i] = 0;

  // We need this for the root grid map
  Enzo_Dims_create(NumberOfRootGrids, TopGridRank, LayoutTemp);
  for (dim = 0; dim < TopGridRank; dim++)
    Layout[TopGridRank-1-dim] = LayoutTemp[dim];
  for (dim = TopGridRank; dim < MAX_DIMENSION; dim++)
    Layout[TopGridRank-1-dim] = 0;

  /* Fill out (processor number)->(grid number) map */

  int gridcount = 0;
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel) {
    Temp->GridData->ReturnGridInfo(&Rank, GridDims, LeftEdge, RightEdge);
    for (dim = 0; dim < Rank; dim++)
      GridPosition[dim] = int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
			      (DomainRightEdge[dim] - DomainLeftEdge[dim]));
    grid_num = GridPosition[0] + 
      Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
    RootProcessors[gridcount] = Temp->GridData->ReturnProcessorNumber();
    GridMap[grid_num] = gridcount++;
  } // ENDFOR level-0 grids

  /* Count level-1 cells contained in each root grid */

  for (Temp = LevelArray[1]; Temp; Temp = Temp->NextGridThisLevel) {
    Temp->GridData->ReturnGridInfo(&Rank, GridDims, LeftEdge, RightEdge);
    for (dim = 0, size = 1; dim < Rank; dim++) {
      GridPosition[dim] = int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
			      (DomainRightEdge[dim] - DomainLeftEdge[dim]));
      size *= GridDims[dim];
    }
    grid_num = GridPosition[0] + 
      Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
    RootGridID = GridMap[grid_num];
    NumberOfSubgridCells[RootGridID] += size;
  } // ENDFOR level-1 grids

  /* Now that we have the number of subgrid cells in each topgrid, we
     can distribute the topgrid tiles so each compute node have a
     similar total number of subgrid cells.  Place them cyclical and
     then move them around. */

  int Node, NewProc;
  int NumberOfNodes = NumberOfProcessors / CoresPerNode;
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
      if (CellsPerNode[i] > MaxVal) {
       MaxVal = CellsPerNode[i];
       MaxNode = i;
      }
      if (CellsPerNode[i] < MinVal) {
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

  delete [] GridMap;
  delete [] NumberOfSubgridCells;
  delete [] CellsPerNode;
  delete [] GridsPerNode;
  delete [] RootProcessors;

  return SUCCESS;
}
