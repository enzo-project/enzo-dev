/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  John Wise (July, 2009): Mode 2/3 -- load balance only
/              within a node.
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdio.h>
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
 
// Function prototypes
 
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void fpcol(float *x, int n, int m, FILE *fptr);
double ReturnWallTime(void);
 
#define LOAD_BALANCE_RATIO 1.05
#define NO_SYNC_TIMING
 
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids, int MoveParticles)
{
 
  if (NumberOfProcessors == 1 || NumberOfGrids <= 1)
    return SUCCESS;
 
  /* Initialize */
 
  int i, GridMemory, NumberOfCells, CellsTotal, Particles, GridsMoved, proc;
  float AxialRatio, GridVolume;
  float *ComputeTime = new float[NumberOfGrids];
  float *ProcessorComputeTime = new float[NumberOfProcessors];
  int *NewProcessorNumber = new int[NumberOfGrids];
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
  moving_count ++;
  out_count ++;
#endif

  double tt0, tt1;
#ifdef SYNC_TIMING
  CommunicationBarrier();
#endif
  tt0 = ReturnWallTime();
 
  GridsMoved = 0;
  for (i = 0; i < NumberOfProcessors; i++)
    ProcessorComputeTime[i] = 0;
 
  /* Compute work for each grid. */
 
  for (i = 0; i < NumberOfGrids; i++) {
    proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
    GridHierarchyPointer[i]->GridData->CollectGridInformation
      (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
    //    ComputeTime[i] = GridMemory; // roughly speaking
    ComputeTime[i] = float(NumberOfCells);
    ProcessorComputeTime[proc] += ComputeTime[i];
    NewProcessorNumber[i] = proc;
  }

 // Mode 1: Load balance over all processors.  Mode 2/3: Load balance
 // only within a node.  Assumes scheduling in blocks (2) or
 // round-robin (3).
 int StartProc, EndProc;
 int proc0 = -1, dproc = -1;
 switch (LoadBalancing) {
 case 1:
   StartProc = 0;
   EndProc = NumberOfProcessors;
   break;
 case 2:  // block scheduling (ranger)
   StartProc = CoresPerNode * (MyProcessorNumber / CoresPerNode);
   EndProc = StartProc + CoresPerNode;
   break;
 case 3:  // round-robin scheduling
   dproc = NumberOfProcessors / CoresPerNode;
   proc0 = MyProcessorNumber % dproc;
   break;
 default:
   if (debug)
     printf("Warning: unknown value for LoadBalancing.  Setting to 1.\n");
   StartProc = 0;
   EndProc = NumberOfProcessors;
   LoadBalancing = 1;
   break;
 } // ENDSWITCH LoadBalancing
 
  /* Transfer grids from heavily-loaded processors. */
 
  int Done = FALSE, MinProc = 0, MaxProc = 0;
  while (!Done) {
 
    /* Find min and max */
 
    float MaxVal = 0, MinVal = huge_number;

    MaxProc = -1;
    MinProc = -1;


    if (LoadBalancing == 1 || LoadBalancing == 2) {
      for (i = StartProc; i < EndProc; i++) {
       if (ProcessorComputeTime[i] > MaxVal) {
         MaxVal = ProcessorComputeTime[i];
         MaxProc = i;
       }
      }
      MinVal = MaxVal;//Sometimes min(ProcessorComputeTime) > huge_number, which causes this loop to fail
      for (i = StartProc; i < EndProc; i++) {
       if (ProcessorComputeTime[i] < MinVal) {
         MinVal = ProcessorComputeTime[i];
         MinProc = i;
       }
      }
    }

    else if (LoadBalancing == 3) {
      for (i = proc0; i < NumberOfProcessors; i += dproc) {
       if (ProcessorComputeTime[i] > MaxVal) {
         MaxVal = ProcessorComputeTime[i];
         MaxProc = i;
       }
      }
      MinVal = MaxVal;//Sometimes min(ProcessorComputeTime) > huge_number, which causes this loop to fail
      for (i = proc0; i < NumberOfProcessors; i += dproc) {
       if (ProcessorComputeTime[i] < MinVal) {
         MinVal = ProcessorComputeTime[i];
         MinProc = i;
       }
      }
    }

    if ((MinProc == -1 || MaxProc == -1) && LoadBalancing == 1)
      fprintf(stderr, "TERRIBLE ERROR [P%"ISYM"]: CommunicationLoadBalance unable to find processors.\n", MyProcessorNumber);
    
    /* Mark a grid transfer if the ratio is large enough. */
 
    if (MaxVal > LOAD_BALANCE_RATIO*MinVal) {
 
      /* Find a grid to transfer. */
 
      for (i = 0; i < NumberOfGrids; i++) {
	//proc = GridHierarchyPointer[i]->GridData->ReturnProcessorNumber();
	proc = NewProcessorNumber[i];
//	if (proc == MaxProc)
//	  printf("P%d / N%d :: grid %d - work = %g, min = %g, max = %g\n",
//		 MyProcessorNumber, MyProcessorNumber/CoresPerNode,
//		 i, ComputeTime[i], MinVal, MaxVal);
	if (proc == MaxProc && ComputeTime[i] < 0.5*(MaxVal-MinVal)) {

//	  printf("\t P%d / N%d :: moving grid %d from %d => %d\n",
//		 MyProcessorNumber, MyProcessorNumber/CoresPerNode,
//		 i, MaxProc, MinProc);
 
	  NewProcessorNumber[i] = MinProc;
	  GridsMoved++;
 
	  /* Update processor compute times. */
 
	  ProcessorComputeTime[MaxProc] -= ComputeTime[i];
	  ProcessorComputeTime[MinProc] += ComputeTime[i];
 
	  break;
	}
      }
 
      /* If we didn't find an appropriate transfer then quit. */
 
      if (i == NumberOfGrids) {
	Done = TRUE;
#ifdef MPI_INSTRUMENTATION
	if (MinVal == 0)
	  timer[3] = MaxVal;
	else
	  timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
      }
    }
    else {
      Done = TRUE;
#ifdef MPI_INSTRUMENTATION
      counter[3]++;
      if (MinVal == 0)
	timer[3] = MaxVal;
      else
	timer[3] = MaxVal/MinVal;
#endif /* MPI_INSTRUMENTATION */
    }
  }
 
#ifdef MPI_INSTRUMENTATION
  moving_pct += float(GridsMoved)/NumberOfGrids;
#endif /* MPI_INSTRUMENTATION */

  /* If mode 2 or 3, gather all moves from other nodes.  Processor
     numbers are only valid on the node.  Zero out other values and
     collect from remote nodes. */

  for (i = 0; i < NumberOfProcessors; i++)
    if (LoadBalancing == 2) {
      if (i < StartProc || i >= EndProc)
       ProcessorComputeTime[i] = 0;
    }
    else if (LoadBalancing == 3) {
      if (i % dproc != proc0)
       ProcessorComputeTime[i] = 0;
    }


  if (LoadBalancing == 2 || LoadBalancing == 3) {
    for (i = 0; i < NumberOfGrids; i++) {

      if (LoadBalancing == 2)
       if (NewProcessorNumber[i] < StartProc || NewProcessorNumber[i] >= EndProc)
         NewProcessorNumber[i] = -1;

      // Trash processor number for remote grids.
      // dproc == # of nodes, proc0 == local node #
      if (LoadBalancing == 3)
       if (NewProcessorNumber[i] % dproc != proc0)
         NewProcessorNumber[i] = -1;

    } // ENDFOR

#ifdef USE_MPI
    CommunicationAllReduceValues(NewProcessorNumber, NumberOfGrids, MPI_MAX);
#endif /* USE_MPI */

  } // ENDIF LoadBalancing == 2 || 3


  /* Now we know where the grids are going, transfer them. */

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (i = 0; i < NumberOfGrids; i++) 
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() !=
	NewProcessorNumber[i])
      GridHierarchyPointer[i]->GridData->
	CommunicationMoveGrid(NewProcessorNumber[i], MoveParticles);

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
  if (MyProcessorNumber == ROOT_PROCESSOR && GridsMoved > 0) {
    tt1 = ReturnWallTime();
    printf("LoadBalance: Number of grids moved = %"ISYM" out of %"ISYM" "
	   "(%lg seconds elapsed)\n", GridsMoved, NumberOfGrids, tt1-tt0);
  }
#ifdef UNUSED
  CommunicationSumValues(ProcessorComputeTime, NumberOfProcessors);
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    printf("LoadBalance (grids=%"ISYM"): \n", NumberOfGrids);
    float norm = ProcessorComputeTime[0];
    for (i = 1; i < NumberOfProcessors; i++)
      norm = max(norm, ProcessorComputeTime[i]);
    for (i = 0; i < NumberOfProcessors; i++)
      ProcessorComputeTime[i] /= max(norm, 1.0e-10);
    // WriteListOfFloats(stdout, NumberOfProcessors, ProcessorComputeTime);
    fpcol(ProcessorComputeTime, NumberOfProcessors, 16, stdout);
  }
#endif 

  delete [] NewProcessorNumber;
  delete [] ComputeTime;
  delete [] ProcessorComputeTime;
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[2] += endtime - starttime;
  counter[2] ++;
#endif /* MPI_INSTRUMENTATION */
 
  return SUCCESS;
}
