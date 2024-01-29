/***********************************************************************
/
/  COMMUNICATION ROUTINE: INITIALIZE
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#include <stdlib.h>
#endif /* USE_MPI */
#ifdef _OPENMP
#include "omp.h"
#endif /* _OPENMP */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sched.h>

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
 
/* function prototypes */
void my_exit(int exit_status);

#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...);
#endif


int CommunicationInitialize(Eint32 *argc, char **argv[])
{
 
#ifdef USE_MPI
 
  /* Initialize MPI and get info. */

  MPI_Arg mpi_rank;
  MPI_Arg mpi_size;
  MPI_Comm comm = MPI_COMM_WORLD;

#ifdef _OPENMP
  MPI_Arg desired = MPI_THREAD_FUNNELED;
  MPI_Arg provided;
  MPI_Arg thread;
  MPI_Init_thread(argc, argv, desired, &provided);
  if (desired != provided) {
    thread = omp_get_thread_num();
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == ROOT_PROCESSOR && thread == 0) {
      fprintf(stderr, "desired = %d, provided = %d\n", desired, provided);
      fprintf(stderr, "WARNING: Cannot get proper OpenMP/MPI setting MPI_THREAD_FUNNELED!\n"
	      "--> Hybrid MPI/OpenMPI mode may fail.\n"
	      "--> Set environment variable MPICH_MAX_THREAD_SAFETY to funneled.\n");
    }
  }
#else
   MPI_Init(argc, argv);
#endif
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_create_errhandler(CommunicationErrorHandlerFn, &CommunicationErrorHandler);
  MPI_Comm_set_errhandler(comm, CommunicationErrorHandler);

  MyProcessorNumber = mpi_rank;
  NumberOfProcessors = mpi_size;
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("MPI_Init: NumberOfProcessors = %"ISYM"\n", NumberOfProcessors);
 
#else /* USE_MPI */
 
  MyProcessorNumber  = 0;
  NumberOfProcessors = 1;
 
#endif /* USE_MPI */

#ifdef _OPENMP
  int CoresPerProcessor;
  int tid, cid;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#pragma omp parallel default(shared) private(tid, cid)
{
  CoresPerProcessor = omp_get_num_threads();
  tid = omp_get_thread_num();
  cid = sched_getcpu();
  NumberOfCores = CoresPerProcessor * NumberOfProcessors;
  printf("Hello from rank: %d, thread: %d, affinity: (core = %d) \n", mpi_rank, tid, cid);
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("MPI_Init: NumberOfCores = %"ISYM"\n", NumberOfCores);
}
#else
  NumberOfCores = NumberOfProcessors;
#endif  
 
  CommunicationTime = 0;
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

#ifndef __APPLE__
#ifdef _OPENMP
  int rank, _thread;
  cpu_set_t coremask;
  char clbuf[7 * CPU_SETSIZE], hnbuf[64];

  memset(clbuf, 0, sizeof(clbuf));
  memset(hnbuf, 0, sizeof(hnbuf));
  (void)gethostname(hnbuf, sizeof(hnbuf));
#pragma omp parallel private(_thread, coremask, clbuf)
  {
    _thread = omp_get_thread_num();
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
//#pragma omp barrier
    //printf("Hello from rank %d, thread %d, on %s. (core affinity = %s)\n",
	  // mpi_rank, _thread, hnbuf, clbuf);
  }
#endif /* _OPENMP */
#endif /* __APPLE__ */
 
  return SUCCESS;
}
 
#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...)
{
  char error_string[1024];
  MPI_Arg length, error_class;
  if (*err != MPI_ERR_OTHER) {
      MPI_Error_class(*err, &error_class);
      MPI_Error_string(error_class, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      MPI_Error_string(*err, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      ENZO_FAIL("MPI communication error.");
  } // ENDIF MPI_ERROR
  return;
}
#endif /* USE_MPI */
 
int CommunicationFinalize()
{
 
#ifdef USE_MPI
  MPI_Errhandler_free(&CommunicationErrorHandler);
  MPI_Finalize();
#endif /* USE_MPI */
 
  return SUCCESS;
}

void CommunicationAbort(int status)
{

#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD,status);
#else
  //  my_exit(status);
#endif

  return;
}
