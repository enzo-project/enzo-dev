/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER STARS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified:   John Wise
/  date:       July, 2009 (modified CTP)
/
/  PURPOSE:
/
************************************************************************/
#ifndef OPTIMIZED_CTP 

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
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
void my_exit(int status);
 
// function prototypes
 
#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_StarMoveList;
#endif
 
 
 
 
int CommunicationTransferStars(grid *GridPointer[], int NumberOfGrids)
{

  if (NumberOfGrids == 1)
    return SUCCESS;
 
  struct StarMoveList {
    int FromGrid;
    int ToGrid[6];
    int NumberToMove[6];
    StarBuffer *Pointer[6];
  };
 
  /* This is a loop that keeps going until no stars were moved. */
 
  int NumberOfStarsMoved, Done = FALSE;
 
  while (!Done) {
 
  NumberOfStarsMoved = 0;
 
  int i, j, grid, GridsToSend = 0;
  StarMoveList *SendList = new StarMoveList[NumberOfGrids];
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    for (i = 0; i < 6; i++) {
      SendList[grid].NumberToMove[i] = 0;
      SendList[grid].Pointer[i]      = NULL;
      SendList[grid].ToGrid[i]       = -1;
    }

  /*
  if ( MyProcessorNumber == 0 ) 
  for (j = 0; j < NumberOfGrids; j++)
    fprintf(stderr, "Pa(%"ISYM") CTP grid[%"ISYM"] = %"ISYM"\n",
                    MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfStars());
  */
 
  /* Generate the list of star moves. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->ReturnProcessorNumber() == MyProcessorNumber) {
      SendList[GridsToSend].FromGrid = grid;
      GridPointer[grid]->CommunicationTransferStars
	(GridPointer, NumberOfGrids, 
	 SendList[GridsToSend].ToGrid, 
	 SendList[GridsToSend].NumberToMove, 
	 SendList[GridsToSend].Pointer, COPY_OUT);
      GridsToSend++;
    }
 
  /* Share the star moves. */
 
  /* Allocate the array to receive subgrids. */
 
  StarMoveList *SharedList = NULL;
 
#ifdef USE_MPI
  int NumberOfSharedGrids = 0;
#endif /* USE_MPI */

#ifdef USE_MPI
  MPI_Status Status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Datatype DataTypeByte = MPI_BYTE;
  MPI_Arg Count;
  MPI_Arg SendCount;
  MPI_Arg RecvCount;
  MPI_Arg Source;
  MPI_Arg Dest;
  MPI_Arg Here;
  MPI_Arg Usize;
  MPI_Arg Xsize;
  MPI_Arg Tag;
  MPI_Arg stat;
#endif /* USE_MPI */
 
  if (NumberOfProcessors > 1) {
 
#ifdef USE_MPI
 
    /* Generate a new MPI type corresponding to the StarMoveList struct. */
 
    if (FirstTimeCalled) {
      Count = sizeof(StarMoveList);
      //  fprintf(stderr, "Size of StarMoveList %"ISYM"\n", Count);
      stat = MPI_Type_contiguous(Count, DataTypeByte, &MPI_StarMoveList);
        if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
      stat = MPI_Type_commit(&MPI_StarMoveList);
        if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
      FirstTimeCalled = FALSE;
    }

    int *SendListCount = new int[NumberOfProcessors];

    MPI_Arg *MPI_SendListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_SendListDisplacements = new MPI_Arg[NumberOfProcessors];

    int jjj;

    for ( jjj=0; jjj<NumberOfProcessors; jjj++) {
      SendListCount[jjj] = 0;
      MPI_SendListCount[jjj] = 0;
      MPI_SendListDisplacements[jjj] = 0;
    }
 
    /* Get counts from each processor to allocate buffers. */
 
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */

    SendCount = 1;
    RecvCount = 1;

    //fprintf(stderr, "GridsToSend %"ISYM" on %"ISYM"\n", GridsToSend, MyProcessorNumber);

    stat = MPI_Allgather(&GridsToSend, SendCount, DataTypeInt, SendListCount, RecvCount, DataTypeInt, MPI_COMM_WORLD);
      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}

/*
    if ( MyProcessorNumber == 0 )
    for ( jjj=0; jjj<NumberOfProcessors; jjj++) {
      fprintf(stderr, "SendListCount %"ISYM" on %"ISYM" \n", SendListCount[jjj], jjj);
    } 
*/

    /* Allocate buffers and generated displacement list. */
 
    for (i = 0; i < NumberOfProcessors; i++) {
      MPI_SendListDisplacements[i] = NumberOfSharedGrids;
      NumberOfSharedGrids += SendListCount[i];
      MPI_SendListCount[i] = SendListCount[i];
    }
    if (NumberOfSharedGrids != NumberOfGrids) {
      ENZO_FAIL("CTP error\n");
    }

/*
    if ( MyProcessorNumber == 0 )
    for ( jjj=0; jjj<NumberOfProcessors; jjj++) {
      fprintf(stderr, "SendListCount %"ISYM" SendListDisp %"ISYM" on %"ISYM" \n", MPI_SendListCount[jjj], MPI_SendListDisplacements[jjj], jjj);
    }
*/


    SharedList = new StarMoveList[NumberOfSharedGrids];
 
    /* Perform sharing operation. */

    Count = GridsToSend;
     
    stat = MPI_Allgatherv(SendList, Count, MPI_StarMoveList, SharedList,
			  MPI_SendListCount, MPI_SendListDisplacements,
			  MPI_StarMoveList, MPI_COMM_WORLD);
      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[9] += endtime-starttime;
    counter[9] ++;
    timer[10] += double(NumberOfSharedGrids);
    GlobalCommunication += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
    delete [] SendListCount;
//    delete [] SendListDisplacements;
    delete [] MPI_SendListCount;
    delete [] MPI_SendListDisplacements;
 
#endif /* USE_MPI */
 
    /* Move stars. */
 
    for (j = 0; j < NumberOfGrids; j++) {
 
      int FromGrid = SharedList[j].FromGrid;
 
      /* Loop over transfer directions. */
 
      for (i = 0; i < 6; i++) {
 
	int ToGrid = max(SharedList[j].ToGrid[i], 0);
	int FromProcessor = GridPointer[FromGrid]->ReturnProcessorNumber();
	int ToProcessor = GridPointer[ToGrid]->ReturnProcessorNumber();
	int TransferSize = SharedList[j].NumberToMove[i];

	NumberOfStarsMoved += SharedList[j].NumberToMove[i];

	if (TransferSize > 0 && FromProcessor != ToProcessor) {
 
        //  fprintf(stderr, "P%"ISYM" I=%"ISYM" XFER SIZE %"ISYM"\n", MyProcessorNumber, i, TransferSize);

#ifdef USE_MPI	
 
	  /* Send stars (Note: this is not good -- the data transfer
	     will not work across heterogeneous machines). */
 
#ifdef MPI_INSTRUMENTATION
          starttime = MPI_Wtime();
#endif  /* MPI_INSTRUMENTATION */

	  Usize = sizeof(StarBuffer);
	  Xsize = TransferSize;

          // fprintf(stderr, "sizeof SendCount %"ISYM"\n", sizeof(SendCount));
          // fprintf(stderr, "sizeof MPI_Arg %"ISYM"\n", sizeof(MPI_Arg));
 
	  if (FromProcessor == MyProcessorNumber) {

//          SendCount = TransferSize*sizeof(StarBuffer);
	    SendCount = Usize * Xsize;
            Dest = ToProcessor;
	    Here = MyProcessorNumber;

	    if( SendCount > 0 ) {
	    // fprintf(stderr, "P%"ISYM" Sending %"ISYM" bytes to %"ISYM" [%"ISYM" x %"ISYM"]\n", Here, SendCount, Dest, Usize, Xsize); 
	    stat = MPI_Send(SharedList[j].Pointer[i], SendCount, DataTypeByte, Dest, MPI_TRANSFERPARTICLE_TAG, MPI_COMM_WORLD);
	      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
	    }
	  }
 
	  /* Receive stars. */
 
	  if (ToProcessor == MyProcessorNumber) {
	    SharedList[j].Pointer[i] = new StarBuffer[TransferSize];

//          RecvCount = TransferSize*sizeof(StarBuffer);
	    RecvCount = Xsize*Usize;
            Source = FromProcessor;
	    Here = MyProcessorNumber;

	    // fprintf(stderr, "P%"ISYM" Post receive %"ISYM" bytes from %"ISYM" [%"ISYM" x %"ISYM"]\n", Here, RecvCount, Source, Usize, Xsize);
	    stat = MPI_Recv(SharedList[j].Pointer[i], RecvCount, DataTypeByte, Source, MPI_TRANSFERPARTICLE_TAG, MPI_COMM_WORLD, &Status);
	      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
	    stat = MPI_Get_count(&Status, DataTypeByte, &RecvCount);
	      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
	    // fprintf(stderr, "P%"ISYM" Received %"ISYM" bytes from %"ISYM" [%"ISYM"]\n", Here, RecvCount, Source, Usize);
	  }
 
#ifdef MPI_INSTRUMENTATION
          double endtime = MPI_Wtime();
	  RecvComm += endtime-starttime;
	  CommunicationTime += endtime-starttime;
#endif  /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
	} // end: if (TransferSize > 0...
 
	/* If this is not on my processor, then clean up pointer, since it
	   doesn't make sense on this processor (or if nothing moved). */
 
	if ((MyProcessorNumber != FromProcessor &&
	     MyProcessorNumber != ToProcessor) || TransferSize == 0)
	  SharedList[j].Pointer[i] = NULL;
 
      } // end: loop over i (directions)
    } // end: loop over j (grids)
 
  } else {
    SharedList = SendList;  // if there is only one processor
    for (grid = 0; grid < NumberOfGrids; grid++)
      for (i = 0; i < 6; i++)
	NumberOfStarsMoved += SharedList[grid].NumberToMove[i];
  }
 
  /* Copy stars back to grids. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->ReturnProcessorNumber() == MyProcessorNumber) {
 
      int LocalNumberToMove[6], counter = 0;
      StarBuffer *LocalPointer[6];
      for (i = 0; i < 6; i++)
	LocalNumberToMove[i] = 0;
 
      /* Extract grid lists to put into this grid. */
 
      for (j = 0; j < NumberOfGrids; j++)
	for (i = 0; i < 6; i++)
	  if (SharedList[j].ToGrid[i] == grid) {
	    LocalNumberToMove[counter] = SharedList[j].NumberToMove[i];
	    LocalPointer[counter++] = SharedList[j].Pointer[i];
	  }
 
      /* Copy stars back (SharedList[grid].ToGrid is a dummy, not used). */
 
      GridPointer[grid]->CommunicationTransferStars
	(GridPointer, NumberOfGrids, SharedList[grid].ToGrid, 
	 LocalNumberToMove, LocalPointer, COPY_IN);
 
    } // end: if grid is on my processor

  /* Set number of particles so everybody agrees. */

#ifdef UNUSED 
  if (NumberOfProcessors > 1) {
    int *Changes = new int[NumberOfGrids];
    for (j = 0; j < NumberOfGrids; j++)
      Changes[j] = 0;
    for (j = 0; j < NumberOfGrids; j++)
      for (i = 0; i < 6; i++)
	if (SharedList[j].ToGrid[i] != -1) {
	  Changes[SharedList[j].FromGrid] -= SharedList[j].NumberToMove[i];
	  Changes[SharedList[j].ToGrid[i]] += SharedList[j].NumberToMove[i];
	}
    for (j = 0; j < NumberOfGrids; j++) {
      if (GridPointer[j]->ReturnProcessorNumber() != MyProcessorNumber)
	GridPointer[j]->SetNumberOfStars(
		         GridPointer[j]->ReturnNumberOfStars()+Changes[j]);
      //      printf("Pb(%"ISYM") CTP grid[%"ISYM"] = %"ISYM"\n", MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfStars());
    }
    delete [] Changes;
  }
#endif /* UNUSED */
 
  /* Set number of stars so everybody agrees. */
 
  if (NumberOfProcessors > 1) {
    int *Changes = new int[NumberOfGrids];
    for (j = 0; j < NumberOfGrids; j++)
      Changes[j] = 0;
    for (j = 0; j < NumberOfGrids; j++)
      for (i = 0; i < 6; i++)
	if (SharedList[j].ToGrid[i] != -1) {
	  Changes[SharedList[j].FromGrid] -= SharedList[j].NumberToMove[i];
	  Changes[SharedList[j].ToGrid[i]] += SharedList[j].NumberToMove[i];
	}
    for (j = 0; j < NumberOfGrids; j++) {
      if (GridPointer[j]->ReturnProcessorNumber() != MyProcessorNumber)
	GridPointer[j]->SetNumberOfStars(
		         GridPointer[j]->ReturnNumberOfStars()+Changes[j]);
      //      printf("Pb(%"ISYM") CTP grid[%"ISYM"] = %"ISYM"\n", MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfParticles());
    }
    delete [] Changes;
  }

  /* CleanUp. */
 
  for (j = 0; j < NumberOfGrids; j++)
    for (i = 0; i < 6; i++) {
      delete [] SharedList[j].Pointer[i];
    }
 
  if (SendList != SharedList)
    delete [] SendList;
 
  delete [] SharedList;
 
  /* Check for completion. */
 
  if (debug)
    printf("CommunicationTransferStars: moved = %"ISYM"\n",
	   NumberOfStarsMoved);
  if (NumberOfStarsMoved == 0)

    Done = TRUE;
 
  }
 
  return SUCCESS;
}
#endif /* OPTIMIZED_CTP */
