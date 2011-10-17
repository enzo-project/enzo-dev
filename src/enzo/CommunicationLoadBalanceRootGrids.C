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
#define MOVES_PER_LOOP 20

int Enzo_Dims_create(int nnodes, int ndims, int *dims);
void icol(int *x, int n, int m, FILE *log_fptr);
void fpcol(Eflt64 *x, int n, int m, FILE *fptr);
int DetermineNumberOfNodes(void);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int LoadBalanceSimulatedAnnealing(int NumberOfGrids, int NumberOfNodes, 
				  int* &ProcessorNumbers, double* NumberOfCells, 
				  double* NumberOfSubcells);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int CommunicationLoadBalanceRootGrids(LevelHierarchyEntry *LevelArray[], 
				      int TopGridRank, int CycleNumber)
{

  if (NumberOfProcessors == 1 || !(LoadBalancing == 2 || LoadBalancing == 3))
    return SUCCESS;

  if (CycleNumber % LoadBalancingCycleSkip != 0)
    return SUCCESS;

  /* If this is a zoom-in calculation, we shouldn't load balance the
     root grids by nodes but the finest static grid. */

  if (StaticRefineRegionLevel[0] != INT_UNDEFINED)
    return SUCCESS;

  /* Declarations */

  LevelHierarchyEntry *Temp;

  char line[MAX_LINE_LENGTH];
  int i, level, dim, Rank, ThisLevel, dummy, ThisTask;
  int size, grid_num, RootGridID, NumberOfRootGrids;
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION], GridPosition[MAX_DIMENSION];
  int GridDims[MAX_DIMENSION];
  double *NumberOfCells, *NumberOfSubgridCells;
  int *RootProcessors, *GridMap;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  int NumberOfNodes = DetermineNumberOfNodes();

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Layout[dim] = 0;
    GridPosition[dim] = 0;
  }
  
  // Root processor does all of the work, and broadcasts the new
  // processor numbers

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  // Count number of level-0 grids
  NumberOfRootGrids = 0;
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
    NumberOfRootGrids++;

  // Allocate some memory
  GridMap = new int[NumberOfRootGrids];
  RootProcessors = new int[NumberOfRootGrids];
  NumberOfCells = new double[NumberOfRootGrids];
  NumberOfSubgridCells = new double[NumberOfRootGrids];

  for (i = 0; i < NumberOfRootGrids; i++) {
    NumberOfCells[i] = 0;
    NumberOfSubgridCells[i] = 0;
  }

  // We need this for the root grid map
  Enzo_Dims_create(NumberOfRootGrids, TopGridRank, LayoutTemp);
  for (dim = 0; dim < TopGridRank; dim++)
    Layout[TopGridRank-1-dim] = LayoutTemp[dim];
  for (dim = TopGridRank; dim < MAX_DIMENSION; dim++)
    Layout[TopGridRank-1-dim] = 0;

  /* Fill out (processor number)->(grid number) map in level-0, then
     count subgrid cells contained in each root grid */

  int gridcount = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {

      Temp->GridData->ReturnGridInfo(&Rank, GridDims, LeftEdge, RightEdge);
      for (dim = 0, size = 1; dim < Rank; dim++) {
	GridPosition[dim] = int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
				(DomainRightEdge[dim] - DomainLeftEdge[dim]));
	size *= GridDims[dim];
      }
      grid_num = GridPosition[0] + 
	Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);

      if (level == 0) {
	RootProcessors[gridcount] = Temp->GridData->ReturnProcessorNumber();
	NumberOfCells[gridcount] += size;
	GridMap[grid_num] = gridcount++;
      } 

      else {
	RootGridID = GridMap[grid_num];
	NumberOfSubgridCells[RootGridID] += size;
      }

    } // ENDFOR grids
  } // ENDFOR levels

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

  } // ENDIF ROOT_PROCESSOR

#ifdef USE_MPI
  MPI_Bcast(&NumberOfRootGrids, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
  if (MyProcessorNumber != ROOT_PROCESSOR)
    RootProcessors = new int[NumberOfRootGrids];
  MPI_Bcast(RootProcessors, NumberOfRootGrids, IntDataType, ROOT_PROCESSOR, 
	    MPI_COMM_WORLD);
#endif

  /* Move the grids to their new processors */

  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, 0, &Grids);

  bool Done = false;
  int nmove, igrid = 0, StartGrid = 0, EndGrid = 0;
  while (!Done) {

    nmove = 0;
    StartGrid = EndGrid;
    while (nmove < MOVES_PER_LOOP && igrid < NumberOfGrids) {
      if (Grids[igrid]->GridData->ReturnProcessorNumber() != RootProcessors[igrid])
	nmove++;
      igrid++;
    }
    EndGrid = igrid;
    if (igrid == NumberOfGrids)
      Done = true;

//    if (debug) 
//      printf("LBRoot: done = %d, grids %d -> %d, nmove = %d\n", 
//	     Done, StartGrid, EndGrid, nmove);

    /* Post receives */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;

    for (i = StartGrid; i < EndGrid; i++)
      if (Grids[i]->GridData->ReturnProcessorNumber() != RootProcessors[i])
	Grids[i]->GridData->CommunicationMoveGrid(RootProcessors[i], TRUE);

    /* Send grids */

    CommunicationDirection = COMMUNICATION_SEND;

    for (i = StartGrid; i < EndGrid; i++)
      if (Grids[i]->GridData->ReturnProcessorNumber() != RootProcessors[i]) {
	if (RandomForcing)  //AK
	  Grids[i]->GridData->AppendForcingToBaryonFields();
	Grids[i]->GridData->CommunicationMoveGrid(RootProcessors[i], TRUE);
      }

    /* Receive grids */

    if (CommunicationReceiveHandler() == FAIL)
      ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
    
    /* Update processor numbers */
    
    for (i = StartGrid; i < EndGrid; i++) {
      Grids[i]->GridData->SetProcessorNumber(RootProcessors[i]);
      if (RandomForcing)  //AK

	Grids[i]->GridData->RemoveForcingFromBaryonFields();
    }

  } // ENDWHILE !Done

  delete [] Grids;
  delete [] RootProcessors;

  return SUCCESS;
}
