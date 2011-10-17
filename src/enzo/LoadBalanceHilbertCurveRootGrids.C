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

#define FUZZY_BOUNDARY 0.1
#define FUZZY_ITERATIONS 10

int LoadBalanceHilbertCurveRootGrids(FLOAT *GridCenters[], int *CellCount,
				     int NumberOfGrids, int* &RootProcessors)
{

  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;

  /* Initialize */
  
  RootProcessors = new int[NumberOfGrids];
  int *GridWork = new int[NumberOfGrids];
  hilbert_data *HilbertData = new hilbert_data[NumberOfGrids];
  int *BlockDivisions = new int[NumberOfProcessors];
  int *ProcessorWork = new int[NumberOfProcessors];

  FLOAT ThisCenter[MAX_DIMENSION];
  int TotalWork, WorkThisProcessor, WorkPerProcessor, WorkLeft;
  int i, dim, grid_num, Rank, block_num, Dims[MAX_DIMENSION];
  int iter;

  /* Compute the position of each grid on a Hilbert curve */
  // TODO: PARALLELIZE

  for (i = 0; i < NumberOfGrids; i++) {

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ThisCenter[dim] = GridCenters[dim][i];

    HilbertData[i].grid_num = i;
    HilbertData[i].hkey = HilbertCurve3D(ThisCenter);
    
  } // ENDFOR grids

  /* Sort the grids along the curve and partition it into pieces with
     equal amounts of work. */

  //qsort(HilbertData, NumberOfGrids, sizeof(hilbert_data), compare_hkey);
  std::sort(HilbertData, HilbertData+NumberOfGrids, cmp_hkey());
  TotalWork = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    GridWork[i] = CellCount[HilbertData[i].grid_num];
    TotalWork += CellCount[HilbertData[i].grid_num];
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

  /* Mark the new processor numbers with the above divisions. */

  block_num = 0;
  for (i = 0; i < NumberOfGrids; i++) {
    if (i > BlockDivisions[block_num]) 
      block_num++;
    RootProcessors[HilbertData[i].grid_num] = block_num;
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
	if (2*GridWork[grid_num] < WorkDifference) {
//	  if (debug)
//	    printf("Moving grid %d (work=%d) from P%d -> P%d\n",
//		   grid_num, GridWork[grid_num], LoadedBlock, UnloadedBlock);
	  ProcessorWork[LoadedBlock] -= GridWork[grid_num];
	  ProcessorWork[UnloadedBlock] += GridWork[grid_num];
	  RootProcessors[HilbertData[grid_num].grid_num] = UnloadedBlock;
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

  /* Cleanup */
  
  delete [] GridWork;
  delete [] HilbertData;
  delete [] ProcessorWork;
  delete [] BlockDivisions;

  return SUCCESS;

}
