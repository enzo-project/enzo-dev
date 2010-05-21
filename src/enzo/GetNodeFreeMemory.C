#ifdef ENABLE_TASKMAP
#ifdef USE_MPI

/***********************************************************************
/
/  DETERMINE NODE FREE MEMORY
/
/  written by: Robert Harkness
/  date:       February, 2006
/
/  PURPOSE:
/
************************************************************************/

#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
void my_exit(int status);


extern "C" Eint64 FreeRealMem(char *name);




int GetNodeFreeMemory(void)
{

  int i, j, k, l, m, n;
  int err;

  Eint64 free;

  int *actual_node;
  int *node_input;
  Eint64 *free_input;
  Eint64 *actual_free;
  double *freemem;

  char node_name[80];


  MPI_Arg nt;
  MPI_Arg id;
  MPI_Arg nn;
  MPI_Arg lname;
  MPI_Arg node_number;
  MPI_Datatype DataTypeInt = (sizeof(Eint64) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  err = MPI_Comm_size(MPI_COMM_WORLD, &nt);
    if( err != 0 ){my_exit(EXIT_FAILURE);}
  err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if( err != 0 ){my_exit(EXIT_FAILURE);}
  err = MPI_Get_processor_name(node_name, &lname);
    if( err != 0 ){my_exit(EXIT_FAILURE);}

  free = FreeRealMem(node_name);

#ifdef DATASTAR
  sscanf(node_name, "ds%3d", &node_number);
#endif

#ifdef TERAGRID
  sscanf(node_name, "tg-c%3d", &node_number);
#endif

#ifdef SEABORG
  sscanf(node_name, "s%5d", &node_number);
#endif

  fprintf(stderr, "Proc %"ISYM" of %"ISYM" is on node %s [%"ISYM"] with %lld MBytes free\n", id, nt, node_name, node_number, free);

  MPI_Barrier(MPI_COMM_WORLD);

  freemem = new double[MAX_NUMBER_OF_NODES];
  actual_node = new int[MAX_NUMBER_OF_NODES];
  node_input = new int[MAX_NUMBER_OF_NODES];
  actual_free = new Eint64[MAX_NUMBER_OF_NODES];
  free_input = new Eint64[MAX_NUMBER_OF_NODES];

  for ( i = 0; i < MAX_NUMBER_OF_NODES; i++ ) {
    NodeMem[i] = 0.0;
    NodeMap[i] = 0;
    actual_node[i] = 0;
    node_input[i] = 0;
    actual_free[i] = 0;
    free_input[i] = 0;
    freemem[i] = 0.0;
  }

//  Get a list of the actual node names to define properties
//  MAX_TASKS_PER_NODE is the number running in this job and NOT
//  necessarily the same as in the target job

  nn = (nt-1)/MAX_TASKS_PER_NODE+1;

  if ( id == 0 )
    fprintf(stderr, "Number of nodes %"ISYM"\n", nn);

  if ( id % MAX_TASKS_PER_NODE == 0 ) {
    node_input[id/MAX_TASKS_PER_NODE] = node_number;
    free_input[id/MAX_TASKS_PER_NODE] = free;
  }

  MPI_Allreduce(node_input, actual_node, nn, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(free_input, actual_free, nn, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

  for ( i = 0; i < nn; i++) {
    NodeMem[i] = ((double) actual_free[i])/1024.0;
    NodeMap[i] = actual_node[i];
  }

  if ( id == 0 ) 
    for ( i = 0; i < nn; i++) {
      fprintf(stderr, "Logical Node %3"ISYM"  NodeMap %5"ISYM"  NodeMem %6.2lf\n", i, NodeMap[i], NodeMem[i]);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  delete [] node_input;
  delete [] free_input;
  delete [] actual_node;
  delete [] actual_free;
  delete [] freemem;

  return SUCCESS;

}
#endif  /* USE_MPI */
#endif  /* ENABLE_TASKMAP */
