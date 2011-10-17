/***********************************************************************
/
/  COMMUNICATION ROUTINE: SUM PARTICLE MASS FLAGGING FIELD
/
/  written by: John Wise
/  date:       May, 2009
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
#define NO_TIMING
 
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
#include "communication.h"
#include "SortCompareFunctions.h"

/* Because we're basically doing a non-blocking MPI_Reduce(MPI_SUM),
   we need to allocate buffers for each on all processors.  Split up
   the work by summing over batches of processors, determined by the
   total number of cells on the level. */

#define CELLS_PER_LOOP 3e8
#define GRIDS_PER_LOOP 100000

void my_exit(int status);
double ReturnWallTime(void);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
int CommunicationBufferPurge(void);
Eint32 compare_proc1(const void *a, const void *b);
Eint32 compare_grid1(const void *a, const void *b);

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_TwoInt;
#endif

int DepositParticleMassFlaggingField(LevelHierarchyEntry* LevelArray[],
				     int level, bool AllLocal)
{

  /* Check if there are any grids on this level, and if we're already
     at the maximum refinement level. */

  if (LevelArray[level] == NULL)
    return SUCCESS;
  
  if (level >= MaximumRefinementLevel)
    return SUCCESS;

  /* Check if there's we're refining by particle mass or must-refine
     particles. */

  int method, ParticleMassMethod;
  bool ParticleFlaggingOn = false;

  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    ParticleFlaggingOn |= (CellFlaggingMethod[method] == 4 ||
			   CellFlaggingMethod[method] == 8);
    if (CellFlaggingMethod[method] == 4)
      ParticleMassMethod = method;
  }

  if (!ParticleFlaggingOn) {
    //fprintf(stdout, "DepositParticleMassFlaggingField not needed!\n"); 
    return SUCCESS;
  }

  LevelHierarchyEntry *Temp;

  if (NumberOfProcessors == 1 || AllLocal) {

    int Zero = 0;

    CommunicationDirection = COMMUNICATION_SEND;
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)

      // Still have to check for the correct processor because this
      // isn't done in the routine because the particles aren't local
      // most of the time.
      if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())
	if (Temp->GridData->SetParticleMassFlaggingField
	    (Zero, Zero, level, ParticleMassMethod) == FAIL) {
	  ENZO_FAIL("Error in grid->SetParticleMassFlaggingField(send).\n");
	}

    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    
  } // ENDIF serial

  else {

#ifdef USE_MPI

    MPI_Arg Count;
    MPI_Arg stat;
    MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
    MPI_Datatype DataTypeByte = MPI_BYTE;
    int sizeof_twoint = sizeof(two_int);

    int *NumberOfSends = new int[NumberOfProcessors];
    int *RecvListCount = new int[NumberOfProcessors];
    MPI_Arg *MPI_SendListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_SendListDisplacements = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_RecvListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_RecvListDisplacements = new MPI_Arg[NumberOfProcessors];

    /* Create a new MPI type corresponding to the structure, two_int */

    if (FirstTimeCalled) {
      Count = sizeof_twoint;
      stat = MPI_Type_contiguous(Count, DataTypeByte, &MPI_TwoInt);
      stat |= MPI_Type_commit(&MPI_TwoInt);
      if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);
      FirstTimeCalled = FALSE;
    }

    HierarchyEntry **Grids;
    int i, count, proc, grid1, StartGrid, EndGrid, NumberOfGrids;
    int TotalNumberOfRecv, TotalNumberOfSends;
    two_int *SendList = NULL;
    two_int *SharedList = NULL;

    /* Make a grid array */

    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

    /* Determine which processors will be sending flagging fields to
       this processor. */

    TotalNumberOfSends = 0;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      if (Grids[grid1]->GridData->ReturnNumberOfParticles() > 0 &&
	  Grids[grid1]->GridData->ReturnProcessorNumber() != MyProcessorNumber)
	TotalNumberOfSends++;

    for (proc = 0; proc < NumberOfProcessors; proc++)
      NumberOfSends[proc] = 0;

    count = 0;
    SendList = new two_int[TotalNumberOfSends];
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      proc = Grids[grid1]->GridData->ReturnProcessorNumber();
      if (Grids[grid1]->GridData->ReturnNumberOfParticles() > 0 &&
	  MyProcessorNumber != proc) {
	SendList[count].grid = grid1;
	SendList[count].proc = proc;
	NumberOfSends[proc]++;
	count++;
      }
    }

    // Fill out displacement lists for MPI_Alltoallv
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      RecvListCount[proc] = 0;
      MPI_RecvListCount[proc] = 0;
      MPI_RecvListDisplacements[proc] = 0;
    }
    
    TotalNumberOfSends = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      MPI_SendListDisplacements[proc] = TotalNumberOfSends;
      TotalNumberOfSends += NumberOfSends[proc];
      MPI_SendListCount[proc] = NumberOfSends[proc];
    }

//    for (i = 0; i < TotalNumberOfSends; i++)
//      printf("P%"ISYM" -- BB -- SendList[%"ISYM"]: %"ISYM" %"ISYM"\n", 
//	     MyProcessorNumber, i, SendList[i].grid, SendList[i].proc);

    // Sort by grid (destination) processor, then replace the
    // processor number with this processor number because the grid
    // processor needs it to post the receives.

#ifdef TIMING
    double t0, t1;
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("--> DPMFF: Before Alltoall, TotalNumberOfSends = %"ISYM"\n",
	     TotalNumberOfSends);
    t0 = ReturnWallTime();
#endif

    //qsort(SendList, TotalNumberOfSends, sizeof_twoint, compare_proc1);
    std::sort(SendList, SendList+TotalNumberOfSends, cmp_proc1());

#ifdef TIMING
    t1 = ReturnWallTime();
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("--> DPMFF: Sorted list by processor in %lg seconds.\n", t1-t0);
    t0 = ReturnWallTime();
#endif

    stat = MPI_Alltoall(NumberOfSends, 1, DataTypeInt,
			RecvListCount, 1, DataTypeInt, MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

    TotalNumberOfRecv = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      MPI_RecvListDisplacements[proc] = TotalNumberOfRecv;
      TotalNumberOfRecv += RecvListCount[proc];
      MPI_RecvListCount[proc] = RecvListCount[proc];
    }

#ifdef TIMING
    t1 = ReturnWallTime();
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("--> DPMFF: After Alltoall, TotalNumberOfRecv = %"ISYM" :: %lg seconds\n",
	     TotalNumberOfRecv, t1-t0);
#endif

    SharedList = new two_int[TotalNumberOfRecv];
    
    /* Perform sharing operation */

    // Since the destination processor really needs the sending
    // processor, replace all of the processor numbers with
    // MyProcessorNumber before sending.

    for (i = 0; i < TotalNumberOfSends; i++)
      SendList[i].proc = MyProcessorNumber;

//    for (i = 0; i < TotalNumberOfSends; i++)
//      printf("P%"ISYM" -- SendList[%"ISYM"]: %"ISYM" %"ISYM"\n", 
//	     MyProcessorNumber, i, SendList[i].grid, SendList[i].proc);

#ifdef TIMING
    t0 = ReturnWallTime();
#endif
    stat = MPI_Alltoallv(SendList, MPI_SendListCount, MPI_SendListDisplacements,
			   MPI_TwoInt,
			 SharedList, MPI_RecvListCount, MPI_RecvListDisplacements,
			   MPI_TwoInt, MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

#ifdef TIMING
    t1 = ReturnWallTime();
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("--> DPMFF: After Alltoallv :: %lg seconds\n", t1-t0);
#endif

//    for (i = 0; i < TotalNumberOfRecv; i++)
//      printf("P%"ISYM" -- SharedList[%"ISYM"]: %"ISYM" %"ISYM"\n", 
//	     MyProcessorNumber, i, SharedList[i].grid, SharedList[i].proc);

    delete [] NumberOfSends;
    delete [] MPI_SendListCount;
    delete [] MPI_SendListDisplacements;

    delete [] RecvListCount;
    delete [] MPI_RecvListCount;
    delete [] MPI_RecvListDisplacements;

    // For easier searching for the sending processors for each grid,
    // sort by grid number
    //qsort(SharedList, TotalNumberOfRecv, sizeof_twoint, compare_grid1);
    std::sort(SharedList, SharedList+TotalNumberOfRecv, cmp_grid1());

    /* Determine how many processors will participate in each
       communication loop to reduce memory usage (see note at the top
       of this file) */

    int dim, size, ProcessorsPerLoop, TotalNumberOfCells;
    int Rank, ThisDims[MAX_DIMENSION];
    FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

    TotalNumberOfCells = 0;
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
      Temp->GridData->ReturnGridInfo(&Rank, ThisDims, Left, Right);
      size = 1;
      for (dim = 0; dim < Rank; dim++)
	size *= ThisDims[dim];
      TotalNumberOfCells += size;
    } // ENDFOR grids
    
    ProcessorsPerLoop = 
      int(CELLS_PER_LOOP / (TotalNumberOfCells/NumberOfProcessors));
    ProcessorsPerLoop = min(max(ProcessorsPerLoop, 1), NumberOfProcessors);
//    ProcessorsPerLoop = NumberOfProcessors;

    int *SendProcs;
    int StartProc, EndProc, nSends;
    double tt0, tt1;

    for (StartProc = 0; StartProc < NumberOfProcessors; 
	 StartProc += ProcessorsPerLoop) {

      count = 0;
      EndProc = min(StartProc + ProcessorsPerLoop, NumberOfProcessors);

      for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {

#ifdef TIMING
	tt0 = ReturnWallTime();
#endif
	EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

	/* Post receives */

	CommunicationReceiveIndex = 0;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	CommunicationDirection = COMMUNICATION_POST_RECEIVE;

	if (TotalNumberOfRecv > 0) {

	  for (grid1 = StartGrid; grid1 < EndGrid && count < TotalNumberOfRecv; 
	       grid1++) {

	    // Create a list of processors from which we're expecting messages.
	    nSends = 0;
	    while (SharedList[count].grid <= grid1) {
	      count++;
	      nSends++;
	      if (count == TotalNumberOfRecv) break;
	    }

	    if (nSends > 0) {
	      SendProcs = new int[nSends];
	      for (i = 0; i < nSends; i++)
		SendProcs[i] = SharedList[count-i-1].proc;
	      if (Grids[grid1]->GridData->SetParticleMassFlaggingField
		  (StartProc, EndProc, level, ParticleMassMethod, SendProcs,
		   nSends) == FAIL) {
		ENZO_FAIL("Error in grid->SetParticleMassFlaggingField"
			"(receive).\n");
	      }
	      delete [] SendProcs;
	    } // ENDIF nSends > 0
	  } // ENDFOR grids
	} // ENDIF TotalNumberOfRecv > 0

	/* Calculate and send field */

	CommunicationDirection = COMMUNICATION_SEND;
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	  if (Grids[grid1]->GridData->SetParticleMassFlaggingField
	      (StartProc, EndProc, level, ParticleMassMethod) == FAIL) {
	    ENZO_FAIL("Error in grid->SetParticleMassFlaggingField(send).\n");
	  }

	/* Wait for the receives and sum field */

	if (CommunicationReceiveHandler() == FAIL)
	  ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
	
#ifdef TIMING
	tt1 = ReturnWallTime();
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  printf("DPMFlag: Finished comm. loop from procs %"ISYM" to %"ISYM", "
		 "grids %"ISYM" to %"ISYM" in %lg seconds.\n", 
		 StartProc, EndProc-1, StartGrid, EndGrid-1, tt1-tt0);
#endif

	CommunicationBufferPurge();

      } // ENDFOR grid batches

    } // ENDFOR processor batches

    // Set communication direction back to its default
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

    delete [] SendList;
    delete [] SharedList;
    delete [] Grids;

#endif /* USE_MPI */

  } // ENDELSE serial

  return SUCCESS;

}

/***********************************************************************/

Eint32 compare_grid1(const void *a, const void *b)
{
  struct two_int *ia = (struct two_int*) a;
  struct two_int *ib = (struct two_int*) b;
  if (ia->grid - ib->grid < 0)
    return -1;
  else if (ia->grid - ib->grid > 0)
    return 1;
  return 0;
}
/***********************************************************************/

Eint32 compare_proc1(const void *a, const void *b)
{
  struct two_int *ia = (struct two_int*) a;
  struct two_int *ib = (struct two_int*) b;
  if (ia->proc - ib->proc < 0)
    return -1;
  else if (ia->proc - ib->proc > 0)

    return 1;
  return 0;
}
