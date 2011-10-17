/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE BY A HILBERT CURVE
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
/  NOTES: For a given level, sort the grids on a 3D Hilbert curve, and
/         then partition the list with equal amounts of work.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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
#include "communication.h"
#include "CommunicationUtilities.h"
#include "SortCompareFunctions.h"

double HilbertCurve3D(FLOAT *coord);
Eint32 compare_hkey(const void *a, const void *b);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
double ReturnWallTime(void);
void fpcol(float *x, int n, int m, FILE *fptr);

#define FUZZY_BOUNDARY 0.1
#define FUZZY_ITERATIONS 10
#define NO_SYNC_TIMING

int LoadBalanceHilbertCurve(HierarchyEntry *GridHierarchyPointer[],
			    int NumberOfGrids, int MoveParticles)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */
  
  int *GridWork = new int[NumberOfGrids];
  int *NewProcessorNumber = new int[NumberOfGrids];
  hilbert_data *HilbertData = new hilbert_data[NumberOfGrids];
  int *BlockDivisions = new int[NumberOfProcessors];
  int *ProcessorWork = new int[NumberOfProcessors];

  int TotalWork, WorkThisProcessor, WorkPerProcessor, WorkLeft;
  int i, dim, grid_num, Rank, block_num, Dims[MAX_DIMENSION];
  FLOAT GridCenter[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT BoundingBox[2][MAX_DIMENSION];
  FLOAT BoundingBoxWidthInv[MAX_DIMENSION];
  float GridVolume, AxialRatio;
  int GridMemory, NumberOfCells, CellsTotal, NumberOfParticles;
  int iter;

  double tt0, tt1;
#ifdef SYNC_TIMING
  CommunicationBarrier();
#endif
  tt0 = ReturnWallTime();

  /* Find the bounding box of the grids */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    BoundingBox[0][dim] = +huge_number;
    BoundingBox[1][dim] = -huge_number;
  }

  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[i]->GridData->
      ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      BoundingBox[0][dim] = min(BoundingBox[0][dim], LeftEdge[dim]);
      BoundingBox[1][dim] = max(BoundingBox[1][dim], RightEdge[dim]);
    }
  } // ENDFOR grids
  
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    BoundingBoxWidthInv[dim] = 1.0/(BoundingBox[1][dim] - BoundingBox[0][dim]);

  /* Compute the position of each grid on a Hilbert curve */
  // TODO: PARALLELIZE

  for (i = 0; i < NumberOfGrids; i++) {

    GridHierarchyPointer[i]->GridData->
      ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);

    // Center relative to the bounding box
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      GridCenter[dim] = 0.5 * (LeftEdge[dim] + RightEdge[dim]);
//      GridCenter[dim] = (0.5 * (LeftEdge[dim] + RightEdge[dim]) - 
//			 BoundingBox[0][dim]) * BoundingBoxWidthInv[dim];

    HilbertData[i].grid_num = i;
    HilbertData[i].hkey = HilbertCurve3D(GridCenter);
    
  } // ENDFOR grids

  /* Sort the grids along the curve and partition it into pieces with
     equal amounts of work. */

  //qsort(HilbertData, NumberOfGrids, sizeof(hilbert_data), compare_hkey);
  std::sort(HilbertData, HilbertData+NumberOfGrids, cmp_hkey());
  TotalWork = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[HilbertData[i].grid_num]->GridData->
      CollectGridInformation(GridMemory, GridVolume, NumberOfCells, 
			     AxialRatio, CellsTotal, NumberOfParticles);
    GridWork[i] = CellsTotal;
    TotalWork += CellsTotal;
  }

  /* Partition into nearly equal workloads */

  grid_num = 0;
  WorkLeft = TotalWork;
  for (i = 0; i < NumberOfProcessors-1; i++) {
    WorkThisProcessor = 0;
    WorkPerProcessor = WorkLeft / (NumberOfProcessors-i);
    while (WorkThisProcessor < WorkPerProcessor) {
      WorkThisProcessor += GridWork[grid_num];
      grid_num++;
    } // ENDWHILE

    // Determine if removing the last grid results in a closer match
    if (grid_num > 1)
      if (ABS(WorkThisProcessor - GridWork[grid_num-1] - WorkPerProcessor) <
	  ABS(WorkThisProcessor - WorkPerProcessor)) {
	WorkThisProcessor -= GridWork[grid_num-1];
	grid_num--;
      }

    // -1 because we advanced grid_num before checking the workload
    BlockDivisions[i] = grid_num-1;  
    ProcessorWork[i] = WorkThisProcessor;
    WorkLeft -= WorkThisProcessor;

  } // ENDFOR processors

  /* Fill in the last entry */

  BlockDivisions[i] = NumberOfGrids-1;
  ProcessorWork[i] = WorkLeft;

//  if (debug) {
//    printf("BlockDivisions = ");
//    for (i = 0; i < NumberOfProcessors; i++)
//      printf("%d ", BlockDivisions[i]);
//    printf("\n");
//  }

  /* Mark the new processor numbers with the above divisions. */

  block_num = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    if (i > BlockDivisions[block_num]) 
      block_num++;
    NewProcessorNumber[HilbertData[i].grid_num] = block_num;
  }

  /* Fuzzy boundaries: There will be some work imbalance if we split
     the grids along the Hilbert curve because the grids where the
     work is split are too large.  So we further refine the balancing
     by moving grids to adjacent processors within some epsilon on the
     Hilbert curve. */

  const double CriticalBalance = 0.1 / NumberOfProcessors;
  double div_hkey, min_hkey, max_hkey, global_min_hkey;
  double hkey_boundary;
  char direction;
  int LoadedBlock, UnloadedBlock, WorkDifference;
  int MinWork, MaxWork;
  float WorkImbalance;

  for (iter = 0; iter < FUZZY_ITERATIONS; iter++) {
    MinWork = 0x7FFFFFFF;
    MaxWork = -1;
    for (i = 0; i < NumberOfProcessors-1; i++) {

      if (BlockDivisions[i] == 0) continue;
      
      /* Hilbert key for the division and boundaries of the curve
	 segment that we will move grids */
      div_hkey = HilbertData[BlockDivisions[i]].hkey;
      if (i == 0)
	min_hkey = div_hkey - FUZZY_BOUNDARY * (div_hkey - HilbertData[0].hkey);
      else
	min_hkey = div_hkey - FUZZY_BOUNDARY * 
	  (div_hkey - HilbertData[BlockDivisions[i-1]].hkey);
      max_hkey = div_hkey + FUZZY_BOUNDARY * 
	(HilbertData[BlockDivisions[i+1]].hkey - div_hkey);
//      if (debug)
//	printf("div%d: Hilbert fuzzy range = %lf -> %lf -> %lf\n",
//	       i, min_hkey, div_hkey, max_hkey);
      
      // Which processor has more work?
      if (ProcessorWork[i] > ProcessorWork[i+1]) {
	LoadedBlock = i;
	UnloadedBlock = i+1;
	hkey_boundary = min_hkey;
	direction = -1;
      } else {
	LoadedBlock = i+1;
	UnloadedBlock = i;
	hkey_boundary = max_hkey;
	direction = +1;
      }

      /* Move grids from the loaded to unloaded processor until we
	 reach the Hilbert key boundary */

      grid_num = BlockDivisions[i] + direction;
      while (direction * (hkey_boundary - HilbertData[grid_num].hkey) > 0) {
	WorkDifference = 
	  ProcessorWork[LoadedBlock] - ProcessorWork[UnloadedBlock];
	if (2*GridWork[grid_num] < WorkDifference &&
	    NewProcessorNumber[HilbertData[grid_num].grid_num] == LoadedBlock) {
//	  if (debug)
//	    printf("Moving grid %d (work=%d) from P%d -> P%d\n",
//		   grid_num, GridWork[grid_num], LoadedBlock, UnloadedBlock);
	  ProcessorWork[LoadedBlock] -= GridWork[grid_num];
	  ProcessorWork[UnloadedBlock] += GridWork[grid_num];
	  NewProcessorNumber[HilbertData[grid_num].grid_num] = UnloadedBlock;
	}
	grid_num += direction;
      } // ENDWHILE move grids
      MinWork = min(MinWork, ProcessorWork[i]);
      MaxWork = max(MaxWork, ProcessorWork[i]);
    } // ENDFOR processors

    MinWork = min(MinWork, ProcessorWork[NumberOfProcessors-1]);
    MaxWork = max(MaxWork, ProcessorWork[NumberOfProcessors-1]);
    WorkImbalance = float(MaxWork - MinWork) / float(MinWork);
    if (WorkImbalance < CriticalBalance)
      break;

  } // ENDFOR iterations

#ifdef UNUSED
  float *ww = new float[NumberOfProcessors];
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("LoadBalance (grids=%"ISYM"): \n", NumberOfGrids);
    float norm = ProcessorWork[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorWork[i]);
    for (i = 0; i < NumberOfProcessors; i++)
      ww[i] = float(ProcessorWork[i]) / max(norm, 1.0e-10);
    // WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
    fpcol(ww, NumberOfProcessors, 16, stdout);
  }
  delete [] ww;
#endif /* UNUSED */

  /* Intermediate cleanup */
  
  delete [] GridWork;
  delete [] HilbertData;
  delete [] ProcessorWork;
  delete [] BlockDivisions;

  /* Now we know where the grids are going, move them! */

  int GridsMoved = 0;

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < NumberOfGrids; i++) 
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[i]) {
      GridHierarchyPointer[i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[i], MoveParticles);
      GridsMoved++;
    }

  /* Send grids */

  CommunicationDirection = COMMUNICATION_SEND;

  for (i = 0; i < NumberOfGrids; i++)
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[i]) {
      if (RandomForcing)  //AK
	GridHierarchyPointer[i]->GridData->AppendForcingToBaryonFields();
      GridHierarchyPointer[i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[i], MoveParticles);
    }

  /* Receive grids */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
  
  /* Update processor numbers */
  
  for (i = 0; i < NumberOfGrids; i++) {
    GridHierarchyPointer[i]->GridData->SetProcessorNumber(NewProcessorNumber[i]);
    if (RandomForcing)  //AK
      GridHierarchyPointer[i]->GridData->RemoveForcingFromBaryonFields();
  }

#ifdef SYNC_TIMING
  CommunicationBarrier();
#endif
  if (debug && GridsMoved > 0) {
    tt1 = ReturnWallTime();
    printf("LoadBalance: Number of grids moved = %"ISYM" out of %"ISYM" "
	   "(%lg seconds elapsed)\n", GridsMoved, NumberOfGrids, tt1-tt0);
  }  

  /* Cleanup */
  
  delete [] NewProcessorNumber;

  return SUCCESS;

}

/************************************************************************
   Version that gives back the processor numbers instead of assigning
   them.  Used in CommunicationPartitionGrid
************************************************************************/

int LoadBalanceHilbertCurve(grid *GridPointers[], int NumberOfGrids, 
			    int* &NewProcessorNumber)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */
  
  NewProcessorNumber = new int[NumberOfGrids];
  int *GridWork = new int[NumberOfGrids];
  hilbert_data *HilbertData = new hilbert_data[NumberOfGrids];
  int *BlockDivisions = new int[NumberOfProcessors];
  int *ProcessorWork = new int[NumberOfProcessors];

  int TotalWork, WorkThisProcessor, WorkPerProcessor, WorkLeft;
  int i, dim, grid_num, Rank, block_num, Dims[MAX_DIMENSION];
  FLOAT GridCenter[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT BoundingBox[2][MAX_DIMENSION];
  FLOAT BoundingBoxWidthInv[MAX_DIMENSION];
  float GridVolume, AxialRatio;
  int GridMemory, NumberOfCells, CellsTotal, NumberOfParticles;
  int iter;

  double tt0, tt1;
  CommunicationBarrier();
  tt0 = ReturnWallTime();

  /* Find the bounding box of the grids */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    BoundingBox[0][dim] = +huge_number;
    BoundingBox[1][dim] = -huge_number;
  }

  for (i = 0; i < NumberOfGrids; i++) {
    GridPointers[i]->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      BoundingBox[0][dim] = min(BoundingBox[0][dim], LeftEdge[dim]);
      BoundingBox[1][dim] = max(BoundingBox[1][dim], RightEdge[dim]);
    }
  } // ENDFOR grids
  
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    BoundingBoxWidthInv[dim] = 1.0/(BoundingBox[1][dim] - BoundingBox[0][dim]);

  /* Compute the position of each grid on a Hilbert curve */
  // TODO: PARALLELIZE

  for (i = 0; i < NumberOfGrids; i++) {

    GridPointers[i]->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);

    // Center relative to the bounding box
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      GridCenter[dim] = 0.5 * (LeftEdge[dim] + RightEdge[dim]);
//      GridCenter[dim] = (0.5 * (LeftEdge[dim] + RightEdge[dim]) - 
//			 BoundingBox[0][dim]) * BoundingBoxWidthInv[dim];

    HilbertData[i].grid_num = i;
    HilbertData[i].hkey = HilbertCurve3D(GridCenter);
    
  } // ENDFOR grids

  /* Sort the grids along the curve and partition it into pieces with
     equal amounts of work. */

  //qsort(HilbertData, NumberOfGrids, sizeof(hilbert_data), compare_hkey);
  std::sort(HilbertData, HilbertData+NumberOfGrids, cmp_hkey());
  TotalWork = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    GridPointers[HilbertData[i].grid_num]->
      CollectGridInformation(GridMemory, GridVolume, NumberOfCells, 
			     AxialRatio, CellsTotal, NumberOfParticles);
    GridWork[i] = CellsTotal;
    TotalWork += CellsTotal;
  }

  /* Partition into nearly equal workloads */

  grid_num = 0;
  WorkLeft = TotalWork;
  for (i = 0; i < NumberOfProcessors-1; i++) {
    WorkThisProcessor = 0;
    WorkPerProcessor = WorkLeft / (NumberOfProcessors-i);
    while (WorkThisProcessor < WorkPerProcessor) {
      WorkThisProcessor += GridWork[grid_num];
      grid_num++;
    } // ENDWHILE

    // Determine if removing the last grid results in a closer match
    if (grid_num > 1)
      if (ABS(WorkThisProcessor - GridWork[grid_num-1] - WorkPerProcessor) <
	  ABS(WorkThisProcessor - WorkPerProcessor)) {
	WorkThisProcessor -= GridWork[grid_num-1];
	grid_num--;
      }

    // -1 because we advanced grid_num before checking the workload
    BlockDivisions[i] = grid_num-1;  
    ProcessorWork[i] = WorkThisProcessor;
    WorkLeft -= WorkThisProcessor;

  } // ENDFOR processors

  /* Fill in the last entry */

  BlockDivisions[i] = NumberOfGrids-1;
  ProcessorWork[i] = WorkLeft;

//  if (debug) {
//    printf("BlockDivisions = ");
//    for (i = 0; i < NumberOfProcessors; i++)
//      printf("%d ", BlockDivisions[i]);
//    printf("\n");
//  }

  /* Mark the new processor numbers with the above divisions. */

  block_num = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    if (i > BlockDivisions[block_num]) 
      block_num++;
    NewProcessorNumber[HilbertData[i].grid_num] = block_num;
  }

  /* Fuzzy boundaries: There will be some work imbalance if we split
     the grids along the Hilbert curve because the grids where the
     work is split are too large.  So we further refine the balancing
     by moving grids to adjacent processors within some epsilon on the
     Hilbert curve. */

  const double CriticalBalance = 0.1 / NumberOfProcessors;
  double div_hkey, min_hkey, max_hkey, global_min_hkey;
  double hkey_boundary;
  char direction;
  int LoadedBlock, UnloadedBlock, WorkDifference;
  int MinWork, MaxWork;
  float WorkImbalance;

  for (iter = 0; iter < FUZZY_ITERATIONS; iter++) {
    MinWork = 0x7FFFFFFF;
    MaxWork = -1;
    for (i = 0; i < NumberOfProcessors-1; i++) {

      if (BlockDivisions[i] == 0) continue;
      
      /* Hilbert key for the division and boundaries of the curve
	 segment that we will move grids */
      div_hkey = HilbertData[BlockDivisions[i]].hkey;
      if (i == 0)
	min_hkey = div_hkey - FUZZY_BOUNDARY * (div_hkey - HilbertData[0].hkey);
      else
	min_hkey = div_hkey - FUZZY_BOUNDARY * 
	  (div_hkey - HilbertData[BlockDivisions[i-1]].hkey);
      max_hkey = div_hkey + FUZZY_BOUNDARY * 
	(HilbertData[BlockDivisions[i+1]].hkey - div_hkey);
//      if (debug)
//	printf("div%d: Hilbert fuzzy range = %lf -> %lf -> %lf\n",
//	       i, min_hkey, div_hkey, max_hkey);
      
      // Which processor has more work?
      if (ProcessorWork[i] > ProcessorWork[i+1]) {
	LoadedBlock = i;
	UnloadedBlock = i+1;
	hkey_boundary = min_hkey;
	direction = -1;
      } else {
	LoadedBlock = i+1;
	UnloadedBlock = i;
	hkey_boundary = max_hkey;
	direction = +1;
      }

      /* Move grids from the loaded to unloaded processor until we
	 reach the Hilbert key boundary */

      grid_num = BlockDivisions[i] + direction;
      while (direction * (hkey_boundary - HilbertData[grid_num].hkey) > 0) {
	WorkDifference = 
	  ProcessorWork[LoadedBlock] - ProcessorWork[UnloadedBlock];
	if (2*GridWork[grid_num] < WorkDifference &&
	    NewProcessorNumber[HilbertData[grid_num].grid_num] == LoadedBlock) {
//	  if (debug)
//	    printf("Moving grid %d (work=%d) from P%d -> P%d\n",
//		   grid_num, GridWork[grid_num], LoadedBlock, UnloadedBlock);
	  ProcessorWork[LoadedBlock] -= GridWork[grid_num];
	  ProcessorWork[UnloadedBlock] += GridWork[grid_num];
	  NewProcessorNumber[HilbertData[grid_num].grid_num] = UnloadedBlock;
	}
	grid_num += direction;
      } // ENDWHILE move grids
      MinWork = min(MinWork, ProcessorWork[i]);
      MaxWork = max(MaxWork, ProcessorWork[i]);
    } // ENDFOR processors

    MinWork = min(MinWork, ProcessorWork[NumberOfProcessors-1]);
    MaxWork = max(MaxWork, ProcessorWork[NumberOfProcessors-1]);
    WorkImbalance = float(MaxWork - MinWork) / float(MinWork);
    if (WorkImbalance < CriticalBalance)
      break;

  } // ENDFOR iterations

#ifdef UNUSED
  float *ww = new float[NumberOfProcessors];
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    printf("LoadBalance (grids=%"ISYM"): \n", NumberOfGrids);
    float norm = ProcessorWork[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorWork[i]);
    for (i = 0; i < NumberOfProcessors; i++)
      ww[i] = float(ProcessorWork[i]) / max(norm, 1.0e-10);
    // WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
    fpcol(ww, NumberOfProcessors, 16, stdout);
  }
  delete [] ww;
#endif /* UNUSED */

  /* Intermediate cleanup */
  
  delete [] GridWork;
  delete [] HilbertData;
  delete [] ProcessorWork;
  delete [] BlockDivisions;

  return SUCCESS;

}

/************************************************************************
   Version where you give the amount of work per grid.  Returns the
   processor numbers instead of assigning them.  Used in
   CommunicationLoadBalancePhotonGrids.
************************************************************************/

int LoadBalanceHilbertCurve(grid *GridPointers[], int NumberOfGrids, 
			    float *InputGridWork, int *NewProcessorNumber)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */
  
  float *GridWork = new float[NumberOfGrids];
  hilbert_data *HilbertData = new hilbert_data[NumberOfGrids];
  int *BlockDivisions = new int[NumberOfProcessors];
  float *ProcessorWork = new float[NumberOfProcessors];

  int GridsLeft;
  int i, dim, grid_num, Rank, block_num, Dims[MAX_DIMENSION];
  FLOAT GridCenter[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT BoundingBox[2][MAX_DIMENSION];
  FLOAT BoundingBoxWidthInv[MAX_DIMENSION];
  float GridVolume, AxialRatio;
  float TotalWork, WorkThisProcessor, WorkPerProcessor, WorkLeft;
  int GridMemory, NumberOfCells, CellsTotal, NumberOfParticles;
  int iter, GridsThisProcessor, GridsPerProcessor;

  double tt0, tt1;
  tt0 = ReturnWallTime();

  /* Find the bounding box of the grids */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    BoundingBox[0][dim] = +huge_number;
    BoundingBox[1][dim] = -huge_number;
  }

  for (i = 0; i < NumberOfGrids; i++) {
    GridPointers[i]->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      BoundingBox[0][dim] = min(BoundingBox[0][dim], LeftEdge[dim]);
      BoundingBox[1][dim] = max(BoundingBox[1][dim], RightEdge[dim]);
    }
  } // ENDFOR grids
  
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    BoundingBoxWidthInv[dim] = 1.0/(BoundingBox[1][dim] - BoundingBox[0][dim]);

  /* Compute the position of each grid on a Hilbert curve */
  // TODO: PARALLELIZE

  for (i = 0; i < NumberOfGrids; i++) {

    GridPointers[i]->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);

    // Center relative to the bounding box
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      GridCenter[dim] = 0.5 * (LeftEdge[dim] + RightEdge[dim]);
//      GridCenter[dim] = (0.5 * (LeftEdge[dim] + RightEdge[dim]) - 
//			 BoundingBox[0][dim]) * BoundingBoxWidthInv[dim];

    HilbertData[i].grid_num = i;
    HilbertData[i].hkey = HilbertCurve3D(GridCenter);
    
  } // ENDFOR grids

  /* Sort the grids along the curve and partition it into pieces with
     equal amounts of work. */

  //qsort(HilbertData, NumberOfGrids, sizeof(hilbert_data), compare_hkey);
  std::sort(HilbertData, HilbertData+NumberOfGrids, cmp_hkey());
  TotalWork = 0;
  for (i = 0; i < NumberOfGrids; i++) {
//    GridPointers[HilbertData[i].grid_num]->
//      CollectGridInformation(GridMemory, GridVolume, NumberOfCells, 
//			     AxialRatio, CellsTotal, NumberOfParticles);
    GridWork[i] = InputGridWork[HilbertData[i].grid_num];
    TotalWork += GridWork[i];
  }

//  if (debug) {
//    printf("TotalWork = %g\n", TotalWork);
//    printf("InputGridWork:\n");
//    fpcol(InputGridWork, NumberOfGrids, 16, stdout);
//  }

  /* Partition into nearly equal workloads */

  grid_num = 0;
  WorkLeft = TotalWork;
  GridsLeft = NumberOfGrids;
  for (i = 0; i < NumberOfProcessors-1; i++) {
    if (WorkLeft == 0) {
      BlockDivisions[i] = grid_num-1;
      ProcessorWork[i] = 0;
      continue;
    }
    WorkThisProcessor = 0;
    GridsThisProcessor = 0;
    WorkPerProcessor = WorkLeft / (NumberOfProcessors-i);
    GridsPerProcessor = GridsLeft / (NumberOfProcessors-i);
    do {
      GridsThisProcessor++;
      GridsLeft--;
      WorkThisProcessor += GridWork[grid_num];
      grid_num++;
    } while (WorkThisProcessor < WorkPerProcessor &&
	     GridsLeft > NumberOfProcessors-i-1);

    // Determine if removing the last grid results in a closer match
    if (GridsThisProcessor > 1)
      if (ABS(WorkThisProcessor - GridWork[grid_num-1] - WorkPerProcessor) <
	  ABS(WorkThisProcessor - WorkPerProcessor)) {
	WorkThisProcessor -= GridWork[grid_num-1];
	grid_num--;
	GridsLeft++;
      }

    // -1 because we advanced grid_num before checking the workload
    BlockDivisions[i] = grid_num-1;  
    ProcessorWork[i] = WorkThisProcessor;
    WorkLeft -= WorkThisProcessor;

  } // ENDFOR processors

  /* Fill in the last entry */

  BlockDivisions[NumberOfProcessors-1] = NumberOfGrids-1;
  ProcessorWork[NumberOfProcessors-1] = WorkLeft;

//  if (debug) {
//    printf("BlockDivisions = ");
//    for (i = 0; i < NumberOfProcessors; i++)
//      printf("%d ", BlockDivisions[i]);
//    printf("\n");
//  }

  /* Mark the new processor numbers with the above divisions. */

  block_num = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    if (i > BlockDivisions[block_num]) 
      block_num++;
    NewProcessorNumber[HilbertData[i].grid_num] = block_num;
  }

  /* Fuzzy boundaries: There will be some work imbalance if we split
     the grids along the Hilbert curve because the grids where the
     work is split are too large.  So we further refine the balancing
     by moving grids to adjacent processors within some epsilon on the
     Hilbert curve. */

  const double CriticalBalance = 0.1 / NumberOfProcessors;
  double div_hkey, min_hkey, max_hkey, global_min_hkey;
  double hkey_boundary;
  char direction;
  int LoadedBlock, UnloadedBlock;
  float WorkDifference, MinWork, MaxWork, WorkImbalance;

#ifdef UNUSED
  for (iter = 0; iter < FUZZY_ITERATIONS; iter++) {
    MinWork = 0x7FFFFFFF;
    MaxWork = -1;
    for (i = 0; i < NumberOfProcessors-1; i++) {

      if (BlockDivisions[i] == 0) continue;
      
      /* Hilbert key for the division and boundaries of the curve
	 segment that we will move grids */
      div_hkey = HilbertData[BlockDivisions[i]].hkey;
      if (i == 0)
	min_hkey = div_hkey - FUZZY_BOUNDARY * (div_hkey - HilbertData[0].hkey);
      else
	min_hkey = div_hkey - FUZZY_BOUNDARY * 
	  (div_hkey - HilbertData[BlockDivisions[i-1]].hkey);
      max_hkey = div_hkey + FUZZY_BOUNDARY * 
	(HilbertData[BlockDivisions[i+1]].hkey - div_hkey);
//      if (debug)
//	printf("div%d: Hilbert fuzzy range = %lf -> %lf -> %lf\n",
//	       i, min_hkey, div_hkey, max_hkey);
      
      // Which processor has more work?
      if (ProcessorWork[i] > ProcessorWork[i+1]) {
	LoadedBlock = i;
	UnloadedBlock = i+1;
	hkey_boundary = min_hkey;
	direction = -1;
      } else {
	LoadedBlock = i+1;
	UnloadedBlock = i;
	hkey_boundary = max_hkey;
	direction = +1;
      }

      /* Move grids from the loaded to unloaded processor until we
	 reach the Hilbert key boundary */

      grid_num = BlockDivisions[i] + direction;
      while (direction * (hkey_boundary - HilbertData[grid_num].hkey) > 0) {
	WorkDifference = 
	  ProcessorWork[LoadedBlock] - ProcessorWork[UnloadedBlock];
	if (2*GridWork[grid_num] < WorkDifference &&
	    NewProcessorNumber[HilbertData[grid_num].grid_num] == LoadedBlock) {
//	  if (debug)
//	    printf("Moving grid %d (work=%d) from P%d -> P%d\n",
//		   grid_num, GridWork[grid_num], LoadedBlock, UnloadedBlock);
	  ProcessorWork[LoadedBlock] -= GridWork[grid_num];
	  ProcessorWork[UnloadedBlock] += GridWork[grid_num];
	  NewProcessorNumber[HilbertData[grid_num].grid_num] = UnloadedBlock;
	}
	grid_num += direction;
      } // ENDWHILE move grids
      MinWork = min(MinWork, ProcessorWork[i]);
      MaxWork = max(MaxWork, ProcessorWork[i]);
    } // ENDFOR processors

    MinWork = min(MinWork, ProcessorWork[NumberOfProcessors-1]);
    MaxWork = max(MaxWork, ProcessorWork[NumberOfProcessors-1]);
    WorkImbalance = float(MaxWork - MinWork) / float(MinWork);
    if (WorkImbalance < CriticalBalance)
      break;

  } // ENDFOR iterations
#endif /* UNUSED */

#ifdef UNUSED
  float *ww = new float[NumberOfProcessors];
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    printf("LoadBalance (grids=%"ISYM"): \n", NumberOfGrids);
    float norm = ProcessorWork[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorWork[i]);
    for (i = 0; i < NumberOfProcessors; i++)
      ww[i] = float(ProcessorWork[i]) / max(norm, 1.0e-10);
    // WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
    fpcol(ww, NumberOfProcessors, 16, stdout);
  }
  delete [] ww;
#endif /* UNUSED */

  /* Intermediate cleanup */
  
  delete [] GridWork;
  delete [] HilbertData;
  delete [] ProcessorWork;
  delete [] BlockDivisions;

  return SUCCESS;

}
