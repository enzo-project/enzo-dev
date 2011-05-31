/***********************************************************************
/
/  COMMUNICATION ROUTINE: DISTRIBUTE STARS IN LIST TO PROCESSORS
/
/  written by: John Wise
/  date:       May, 2009
/  modified:   July, 2009 by John Wise -- adapted for stars
/
/  PURPOSE: Takes a list of star moves and sends/receives stars
/           to all processors
/
************************************************************************/

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
#include "SortCompareFunctions.h"

void my_exit(int status);

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_StarMoveList;
#endif

Eint32 compare_star_proc(const void *a, const void *b);
Eint32 compare_star_grid(const void *a, const void *b);

int CommunicationShareStars(int *NumberToMove, star_data* &SendList,
			    int &NumberOfReceives, star_data* &SharedList)
{

  int i, proc;
  int NumberOfSends;

  // In order to call Alltoallv, we need to sort the list by
  // destination processor

  int TotalNumberToMove = 0;
  int star_data_size = sizeof(star_data);
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];
  //qsort(SendList, TotalNumberToMove, star_data_size, compare_star_proc);
  std::sort(SendList, SendList+TotalNumberToMove, cmp_star_proc());

  SharedList = NULL;

#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Datatype DataTypeByte = MPI_BYTE;
  MPI_Arg Count;
  MPI_Arg SendCount;
  MPI_Arg RecvCount;
  MPI_Arg stat;
#endif /* USE_MPI */

  if (NumberOfProcessors > 1) {

#ifdef USE_MPI

    /* Generate a new MPI type corresponding to the StarMoveList data. */
 
    if (FirstTimeCalled) {
      Count = sizeof(star_data);
      //  fprintf(stderr, "Size of StarMoveList %"ISYM"\n", Count);
      stat = MPI_Type_contiguous(Count, DataTypeByte, &MPI_StarMoveList);
      stat |= MPI_Type_commit(&MPI_StarMoveList);
      if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);
      FirstTimeCalled = FALSE;
    }

    /* Get counts from each processor to allocate buffers. */

    MPI_Arg *MPI_SendListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_SendListDisplacements = new MPI_Arg[NumberOfProcessors];

    int *RecvListCount = new int[NumberOfProcessors];
    MPI_Arg *MPI_RecvListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_RecvListDisplacements = new MPI_Arg[NumberOfProcessors];

    int jjj;

    for (jjj = 0; jjj < NumberOfProcessors; jjj++) {
      RecvListCount[jjj] = 0;
      MPI_RecvListCount[jjj] = 0;
      MPI_RecvListDisplacements[jjj] = 0;
    }

    NumberOfSends = 0;
    for (jjj = 0; jjj < NumberOfProcessors; jjj++) {
      MPI_SendListDisplacements[jjj] = NumberOfSends;
      NumberOfSends += NumberToMove[jjj];
      MPI_SendListCount[jjj] = NumberToMove[jjj];
    }

    SendCount = 1;
    RecvCount = 1;

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */

    /**************************
       Share the star counts
    ***************************/
    
    stat = MPI_Alltoall(NumberToMove, SendCount, DataTypeInt,
			RecvListCount, RecvCount, DataTypeInt, MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

    /* Allocate buffers and generated displacement list. */

    NumberOfReceives = 0;  
    for (i = 0; i < NumberOfProcessors; i++) {
      MPI_RecvListDisplacements[i] = NumberOfReceives;
      NumberOfReceives += RecvListCount[i];
      MPI_RecvListCount[i] = RecvListCount[i];
    }

    SharedList = new star_data[NumberOfReceives];
 
    /******************************
          Share the stars
    ******************************/

    stat = MPI_Alltoallv(SendList, MPI_SendListCount, MPI_SendListDisplacements,
			   MPI_StarMoveList,
			 SharedList, MPI_RecvListCount, MPI_RecvListDisplacements,
			   MPI_StarMoveList,
			 MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);

#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[9] += endtime-starttime;
    counter[9] ++;
    timer[10] += double(NumberOfReceives);
    GlobalCommunication += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
    delete [] MPI_SendListCount;
    delete [] MPI_SendListDisplacements;

    delete [] RecvListCount;
    delete [] MPI_RecvListCount;
    delete [] MPI_RecvListDisplacements;

#endif /* USE_MPI */    

  } // ENDIF multi-processor
  else {
    NumberOfReceives = TotalNumberToMove;
    SharedList = SendList;
  }
  
  // First sort the list by destination grid, so the searching for
  // grids is more efficient.
  //qsort(SharedList, NumberOfReceives, star_data_size, compare_star_grid);
  std::sort(SharedList, SharedList+NumberOfReceives, cmp_star_grid());

  return SUCCESS;

}
