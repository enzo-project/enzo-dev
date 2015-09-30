/***********************************************************************
/
/  COMMUNICATION ROUTINE: DISTRIBUTE PARTICLES IN LIST TO PROCESSORS
/
/  written by: John Wise
/  date:       May, 2009
/  modified:   
/
/  PURPOSE: Takes a list of particle moves and sends/receives particles
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
static MPI_Datatype MPI_ParticleMoveList;
#endif

Eint32 compare_proc(const void *a, const void *b);
Eint32 compare_grid(const void *a, const void *b);

int CommunicationShareParticles(int *NumberToMove, particle_data* &SendList,
				int &NumberOfReceives,
				particle_data* &SharedList)
{

  int i, proc;
  int NumberOfSends;

  // In order to call Alltoallv, we need to sort the list by
  // destination processor

  int TotalNumberToMove = 0;
  int particle_data_size = sizeof(particle_data);
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];
  //qsort(SendList, TotalNumberToMove, particle_data_size, compare_proc);
  std::sort(SendList, SendList+TotalNumberToMove, cmp_proc());

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

    /* Generate a new MPI type corresponding to the ParticleMoveList data. */
 
    if (FirstTimeCalled) {
      Count = sizeof(particle_data);
      //  fprintf(stderr, "Size of ParticleMoveList %"ISYM"\n", Count);
      stat = MPI_Type_contiguous(Count, DataTypeByte, &MPI_ParticleMoveList);
      stat |= MPI_Type_commit(&MPI_ParticleMoveList);
      if (stat != MPI_SUCCESS) ENZO_FAIL("");
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

    /******************************
       Share the particle counts
    ******************************/
    
    stat = MPI_Alltoall(NumberToMove, SendCount, DataTypeInt,
			RecvListCount, RecvCount, DataTypeInt, MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) ENZO_FAIL("");

    /* Allocate buffers and generated displacement list. */

    NumberOfReceives = 0;  
    for (i = 0; i < NumberOfProcessors; i++) {
      MPI_RecvListDisplacements[i] = NumberOfReceives;
      NumberOfReceives += RecvListCount[i];
      MPI_RecvListCount[i] = RecvListCount[i];
    }

    SharedList = new particle_data[NumberOfReceives];
 
    /******************************
          Share the particles
    ******************************/

    stat = MPI_Alltoallv(SendList, MPI_SendListCount, MPI_SendListDisplacements,
			   MPI_ParticleMoveList,
			 SharedList, MPI_RecvListCount, MPI_RecvListDisplacements,
			   MPI_ParticleMoveList,
			 MPI_COMM_WORLD);
    if (stat != MPI_SUCCESS) ENZO_FAIL("");

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
  //qsort(SharedList, NumberOfReceives, particle_data_size, compare_grid);
  std::sort(SharedList, SharedList+NumberOfReceives, cmp_grid());

  return SUCCESS;

}
