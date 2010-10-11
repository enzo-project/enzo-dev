/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS BY RAY TRACING WORK
/
/  written by: John Wise
/  date:       September, 2010
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
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
#include "communication.h"
#include "CommunicationUtilities.h"

#define LOAD_BALANCE_RATIO 1.05
#define MIN_LEVEL 1

/* function prototypes */

int LoadBalanceHilbertCurve(grid *GridPointers[], int NumberOfGrids, 
			    float *InputGridWork, int *NewProcessorNumber);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
double ReturnWallTime(void);
int FindField(int field, int farray[], int numfields);
void fpcol(Eflt32 *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);

Eint32 _compare(const void * a, const void * b)
{
  return ( *(Eflt32*)a - *(Eflt32*)b );
}

int CommunicationLoadBalancePhotonGrids(HierarchyEntry **Grids[], int *NumberOfGrids)
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef TRANSFER // nee TRANSFER

  if (RadiativeTransfer == FALSE)
    return SUCCESS;

  if (RadiativeTransferLoadBalance == FALSE)
    return SUCCESS;

  double tt0, tt1;
#ifdef SYNC_TIMING
  CommunicationBarrier();
#endif
  tt0 = ReturnWallTime();

  /* Initialize */

  float *NonZeroComputeTime, *ComputeTime;
  int *NewProcessorNumber, *NonZeroList;
  grid **NonZeroGrids;
  char *MoveFlag;

  int i, index, index2, lvl, dim, proc, GridsMoved, TotalNumberOfGrids;
  int NumberOfBaryonFields, Nonzero;
  int FieldTypes[MAX_NUMBER_OF_BARYON_FIELDS];
  float AxialRatio, GridVolume;
  float *ProcessorComputeTime = new float[NumberOfProcessors];
  Eflt32 *SortedComputeTime = new Eflt32[NumberOfProcessors];

  GridsMoved = 0;

  TotalNumberOfGrids = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    TotalNumberOfGrids += NumberOfGrids[lvl];

  /* Allocate memory */

  ComputeTime = new float[TotalNumberOfGrids];
  MoveFlag = new char[TotalNumberOfGrids];
  NewProcessorNumber = new int[TotalNumberOfGrids];
  
  /* Compute work for each grid. */

  NumberOfBaryonFields = Grids[0][0]->GridData->
    ReturnNumberOfBaryonFields();
  Grids[0][0]->GridData->ReturnFieldType(FieldTypes);
  int RaySegNum = FindField(RaySegments, FieldTypes, NumberOfBaryonFields);

  index = 0;
  Nonzero = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++) {

      proc = Grids[lvl][i]->GridData->ReturnProcessorNumber();
      Grids[lvl][i]->GridData->SetOriginalProcessorNumber(proc);
      if (MyProcessorNumber == proc) {
	ComputeTime[index] = Grids[lvl][i]->GridData->
	  ReturnTotalNumberOfRaySegments(RaySegNum);
	if (ComputeTime[index] > 0) Nonzero++;
      } else
	ComputeTime[index] = 0.0;

      MoveFlag[index] = FALSE;
      NewProcessorNumber[index] = proc;

  } // ENDFOR grids

  /* Get total compute time over all processors */

  CommunicationSumValues(ComputeTime, TotalNumberOfGrids);
  CommunicationSumValues(&Nonzero, 1);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  /* Find grids with non-zero amount of work and only load balance
     those grids.  Construct a contiguous grid list with work. */

  NonZeroComputeTime = new float[Nonzero];
  NonZeroList = new int[Nonzero];
  for (i = 0, index = 0; i < TotalNumberOfGrids; i++)
    if (ComputeTime[i] > 0) {
      NonZeroComputeTime[index] = ComputeTime[i];
      NonZeroList[index++] = i;
    } // ENDIF

  NonZeroGrids = new grid*[Nonzero];
  index = 0;
  index2 = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++)
      if (ComputeTime[index] > 0)
	NonZeroGrids[index2++] = Grids[lvl][i]->GridData;

  /* Transfer grids from heavily-loaded processors. */

  int *NonZeroProcs = new int[Nonzero];
  LoadBalanceHilbertCurve(NonZeroGrids, Nonzero, NonZeroComputeTime,
			  NonZeroProcs);

  // Change processor number in list.
  for (i = 0; i < Nonzero; i++)
    NewProcessorNumber[NonZeroList[i]] = NonZeroProcs[i];

  delete [] NonZeroComputeTime;
  delete [] NonZeroList;
  delete [] NonZeroGrids;
  delete [] NonZeroProcs;

  } // ENDIF ROOT_PROCESSOR

    for (i = 0; i < NumberOfProcessors; i++) {
      Mean += ProcessorComputeTime[i];
      if (ProcessorComputeTime[i] > MaxVal) {
	MaxVal = ProcessorComputeTime[i];
	MaxProc = i;
      }
    }
    for (i = 0; i < NumberOfProcessors; i++) {
      if (ProcessorComputeTime[i] < MinVal) {
	MinVal = ProcessorComputeTime[i];
	MinProc = i;
      }
    }

    Mean /= NumberOfProcessors;

    for (i = 0; i < NumberOfProcessors; i++)
      SortedComputeTime[i] = ProcessorComputeTime[i];
    qsort(SortedComputeTime, (size_t) NumberOfProcessors, sizeof(Eflt32), 
	  _compare);
    FirstNonZero = NumberOfProcessors-1;
    for (i = 0; i < NumberOfProcessors; i++)
      if (SortedComputeTime[i] > 0) {
	FirstNonZero = i;
	break;
      }
    Median = SortedComputeTime[FirstNonZero + (NumberOfProcessors-FirstNonZero)/2];

    if (MaxVal > LOAD_BALANCE_RATIO*MinVal) {
      /* Find a grid to transfer. */

      for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++) {
      BreakLoop = false;
      for (i = 0; i < NumberOfGrids[lvl]; i++) {
	//proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
	proc = NewProcessorNumber[lvl][i];
	//if (ProcessorComputeTime[proc] > Mean && debug && ComputeTime[lvl][i] > 0)
	//  printf("P%d / L%d :: grid %d - work = %g, min = %g, max = %g, median = %g\n",
	//	 MyProcessorNumber, lvl, i, ComputeTime[lvl][i], MinVal, MaxVal, Median);
	if (ProcessorComputeTime[proc] > Mean && 
	    ComputeTime[lvl][i] < Mean && ComputeTime[lvl][i] > 0 && 
	    MoveFlag[lvl][i] == FALSE) {

//	  if (debug)
//	    printf("\t P%d / L%d :: moving grid %d from %d => %d\n",
//		   MyProcessorNumber, lvl, i, proc, MinProc);
 
	  NewProcessorNumber[lvl][i] = MinProc;
	  GridsMoved++;
 
	  /* Update processor compute times. */
 
	  ProcessorComputeTime[proc] -= ComputeTime[lvl][i];
	  ProcessorComputeTime[MinProc] += ComputeTime[lvl][i];
	  MoveFlag[lvl][i] = TRUE;

	  BreakLoop = true;
	  break;
	}
      } // ENDFOR grids

      if (BreakLoop) break;
 
      /* If we didn't find an appropriate transfer then quit. */
 
      if (i == NumberOfGrids[lvl])
	Done = TRUE;
      
      } // ENDFOR levels

    } // ENDIF !( load balanced )
    else {
      Done = TRUE;
    }
    
  } // ENDWHILE !Done

//  if (debug) {
//    printf("After load balancing:\n");
//    fpcol(ProcessorComputeTime, NumberOfProcessors, 8, stdout);
//  }


  /* Now we know where the grids are going, transfer them. */

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  index = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++) 
      if (Grids[lvl][i]->GridData->ReturnProcessorNumber() !=
	  NewProcessorNumber[index]) {
      Grids[lvl][i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[index], FALSE, FALSE);
      MoveFlag[index] = TRUE;
      GridsMoved++;
    } else {
      MoveFlag[index] = FALSE;
    }

  /* Send grids */

  CommunicationDirection = COMMUNICATION_SEND;

  index = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++)
    if (Grids[lvl][i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[index]) {
      if (RandomForcing)  //AK
	Grids[lvl][i]->GridData->AppendForcingToBaryonFields();
      Grids[lvl][i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[index], FALSE, FALSE);
    }

  /* Receive grids */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

  /* Update processor numbers */

  index = 0;
  for (lvl = MIN_LEVEL; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    for (i = 0; i < NumberOfGrids[lvl]; i++, index++) {
      Grids[lvl][i]->GridData->SetProcessorNumber(NewProcessorNumber[index]);
      if (RandomForcing)  //AK
	Grids[lvl][i]->GridData->RemoveForcingFromBaryonFields();
  }

  CommunicationBarrier();

  if (MyProcessorNumber == ROOT_PROCESSOR && GridsMoved > 0) {
    tt1 = ReturnWallTime();
    printf("PhotonLoadBalance: Number of grids moved = %"ISYM" out of %"ISYM" "
	   "(%lg seconds elapsed)\n", GridsMoved, TotalNumberOfGrids, tt1-tt0);
  }

  /* Cleanup */

  delete [] ComputeTime;
  delete [] MoveFlag;
  delete [] NewProcessorNumber;

#endif /* TRANSFER */

  return SUCCESS;
}
