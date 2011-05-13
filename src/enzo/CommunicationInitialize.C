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

int CommunicationInitialize(Eint32 *argc, char **argv[])
{
 
#ifdef USE_MPI
 
  /* Initialize MPI and get info. */

  MPI_Arg mpi_rank;
  MPI_Arg mpi_size;

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
#pragma omp parallel
  CoresPerProcessor = omp_get_num_threads();
  NumberOfCores = CoresPerProcessor * NumberOfProcessors;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("MPI_Init: NumberOfCores = %"ISYM"\n", NumberOfCores);
#else
  NumberOfCores = NumberOfProcessors;
#endif  
 
  CommunicationTime = 0;
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  return SUCCESS;
}
 
 
 
 
int CommunicationFinalize()
{
 
#ifdef USE_MPI
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
