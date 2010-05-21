#ifdef ENABLE_TASKMAP
#ifdef USE_MPI

/***********************************************************************
/
/  GENERATE GRID TO MPI TASK MAP
/
/  written by: Robert Harkness
/  date:       February, 2006
/  modified:   May 31st, 2007
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




int AssignGridToTaskMap(Eint64 GridIndex[], Eint64 Memory[], int Ntask)
{

  int i, j, k, l, m, n;
  int err;
  int node;
  int gtemp;
  int mnode;
  double maxmem;

  Eint64 free;

  double *grid_size;
  int *grid;
  int *task_node;
  int *task;

  double *freemem;
  int *freecpu;
  int *actual_node;
  int *node_input;
  Eint64 *free_input;
  Eint64 *actual_free;
  int *tasks_per_node;
  double *nodemem;

  char node_name[80];

  grid_size = new double[MAX_NUMBER_OF_TASKS];
  grid = new int[MAX_NUMBER_OF_TASKS];
  task_node = new int[MAX_NUMBER_OF_TASKS];
  task = new int[MAX_NUMBER_OF_TASKS];

  freemem = new double[MAX_NUMBER_OF_NODES];
  freecpu = new int[MAX_NUMBER_OF_NODES];
  actual_node = new int[MAX_NUMBER_OF_NODES];
  node_input = new int[MAX_NUMBER_OF_NODES];
  actual_free = new Eint64[MAX_NUMBER_OF_NODES];
  free_input = new Eint64[MAX_NUMBER_OF_NODES];
  tasks_per_node = new int[MAX_NUMBER_OF_NODES];
  nodemem = new double[MAX_NUMBER_OF_NODES];

  int task_table[MAX_NUMBER_OF_NODES][MAX_TASKS_PER_NODE];
  int task_counter[MAX_NUMBER_OF_NODES];


  for ( i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    grid_size[i] = 0.0;
    task_node[i] = -1;
  }

  for ( n = 0; n < MAX_NUMBER_OF_NODES; n++ ) {
    actual_node[n] = 0;
    actual_free[n] = 0;
    node_input[n] = 0;
    free_input[n] = 0;
    freemem[n] = 0.0;
    freecpu[n] = 0;
    nodemem[n] = 0.0;

    task_counter[n] = 0;
  }

  for ( n = 0; n < MAX_NUMBER_OF_NODES; n++ ) {
    for ( i = 0; i < MAX_TASKS_PER_NODE; i++ ) {
      task_table[n][i] = 0;
    }
  }

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

  fprintf(stderr, "MPI Task %"ISYM" of %"ISYM" is on node %s [%"ISYM"] with %lld MBytes free\n", id, nt, node_name, node_number, free);

  MPI_Barrier(MPI_COMM_WORLD);

//  Get a list of the actual node names to define properties
//  MAX_TASKS_PER_NODE is the number running in this job and NOT
//  necessarily the same as in the target job

  nn = (nt-1)/MAX_TASKS_PER_NODE+1;

  int nnx = nt / MAX_TASKS_PER_NODE;
  int nnr = nt % MAX_TASKS_PER_NODE;

  nn = nnx;
  if (nnr != 0 ) nn = nnx + 1;

  for ( i = nn*MAX_TASKS_PER_NODE; i < (nn+1)*MAX_TASKS_PER_NODE; i++ ) {
    task_node[i] = nn;
  }


  if ( id == 0 )
    fprintf(stderr, "Number of nodes %"ISYM"\n", nn);

  if ( id % MAX_TASKS_PER_NODE == 0 ) {
    node_input[id/MAX_TASKS_PER_NODE] = node_number;
    free_input[id/MAX_TASKS_PER_NODE] = free;
  }

  MPI_Allreduce(node_input, actual_node, nn, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(free_input, actual_free, nn, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

// Define node properties
// As long as the total number of freecpu[] is > nt this should work
// MPI tasks numbers are known for each node
// Make a table of MPI tasks owned by each node

  for ( n = 0; n < nn; n++) {

    tasks_per_node[n] = MAX_TASKS_PER_NODE;
    freecpu[n] = tasks_per_node[n];
    freemem[n] = ((double) actual_free[n])/1024.0;

    for ( j = n*MAX_TASKS_PER_NODE; j < min((n+1)*MAX_TASKS_PER_NODE, nt); j++ ) {
      task_table[n][j-n*MAX_TASKS_PER_NODE] = j;
    }
  }

  if ( id == 0 ) 
    for ( n = 0; n < nn; n++) {
      fprintf(stderr, "Logical Node %3d  Actual Node %5"ISYM"  Initial Memory %6.2lf  Present Memory %6.2lf : ",
                      n, actual_node[n], NodeMem[n], freemem[n]);
      for ( j = 0; j < MAX_TASKS_PER_NODE; j++ ) {
        fprintf(stderr, " %4d", task_table[n][j]);
      }
      fprintf(stderr, "\n");
    }

// Use the initial memory for the distribution!

  for ( n = 0; n < nn; n++) {
    if( actual_node[n] != NodeMap[n] ){my_exit(EXIT_FAILURE);} 
    freemem[n] = NodeMem[n];
    tasks_per_node[n] = MAX_TASKS_PER_NODE;
    freecpu[n] = tasks_per_node[n];
  }


// Generate mapping

  for ( i = 0; i < nt; i++ ) {
    grid[i] = GridIndex[i];
    gtemp = Memory[i];
    grid_size[i] = ((double) gtemp) / 1073741824.0;  // in GBytes
  }

  if ( id == 0 ) {
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "%"ISYM"  %"ISYM"  %6.2f\n", i, grid[i], grid_size[i]);
    }
  }

// Sort grid[] and grid_size[] to monotonically decreasing grid_size

  int last1, save;
  double dsave;

  for ( m = 1; m < nt; m++ ) {
    last1 = nt-m+1;
    for ( l = 1; l < last1; l++) {
      if( grid_size[l] > grid_size[l-1] ) {
        dsave = grid_size[l-1];
        grid_size[l-1] = grid_size[l];
        grid_size[l] = dsave;
        save = grid[l-1];
        grid[l-1] = grid[l];
        grid[l] = save;
      }
    }
  }

  if ( id == 0 ) {
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "%"ISYM"  %"ISYM"  %6.2f\n", i, grid[i], grid_size[i]);
    }
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

/* Here we have on every task
   grid[i]               the grid number
   grid_size[i]          the size of the grid in memory
                         strictly monotonic decreasing

   nt              number of MPI tasks
   nn              number of physical nodes
   freemem[n]      free memory on node n
   freecpu[n]      free cpus on node n
*/

// assign grid[i] to task[j] on node[m]

  for ( i = 0; i < nt; i++ ) {

  maxmem = 0.0;
  mnode = -1;

  for ( n = 0; n < nn; n++ ) {
    if ( freecpu[n] > 0 ) {
      if ( freemem[n] > maxmem ) {
        maxmem = freemem[n];
        mnode = n;
      }
    }
  }

  freecpu[mnode] = freecpu[mnode] - 1;
  freemem[mnode] = freemem[mnode] - grid_size[i];
  task_node[i] = mnode;

  if ( id == 0 )
    if ( freemem[mnode] < 0.0 )
      fprintf(stderr, "memory < 0 task %"ISYM" mnode %"ISYM" freemem %6.2f\n", i, mnode, freemem[mnode]);

  }

  MPI_Barrier(MPI_COMM_WORLD);

  for ( i = 0; i < nt; i++ ) {

     // task[i] = i; // WRONG!

     // Grid 586 assigned to Node 250 is correct
     // which tasks are assigned to Node 250?
     // this is the next task on that node

    task[i] = task_node[i]*MAX_TASKS_PER_NODE + task_counter[task_node[i]];
    task_counter[task_node[i]]++;

    if ( id == 0 ) fprintf(stderr, "Task %"ISYM"  Grid %"ISYM"  Node %"ISYM"  Size %6.2f\n", task[i], grid[i], task_node[i], grid_size[i]);
  }
  if ( id == 0 ) fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");

  MPI_Barrier(MPI_COMM_WORLD);

  for ( m = 1; m < nt; m++ ) {
    last1 = nt-m+1;
    for ( l = 1; l < last1; l++) {
      if( grid[l] < grid[l-1] ) {
        save = grid[l-1];
        grid[l-1] = grid[l];
        grid[l] = save;
        save = task_node[l-1];
        task_node[l-1] = task_node[l];
        task_node[l] = save;
        save = task[l-1];
        task[l-1] = task[l];
        task[l] = save;
        dsave = grid_size[l-1];
        grid_size[l-1] = grid_size[l];
        grid_size[l] = dsave;
      }
    }
  }

  for ( i = 0; i < nt; i++ ) {
    TaskMap[i] = task[i];
    TaskMemory[i] = 0;
  }

  if ( id == 0 ) {
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "Grid %"ISYM"  Task %"ISYM"  Node %"ISYM"  Size %6.2f\n", grid[i], task[i], task_node[i], grid_size[i]);
    }
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for ( i = 0; i < nt; i++ ) {
    nodemem[task_node[i]] = nodemem[task_node[i]] + grid_size[i];
  }

  if ( id == 0 )
  for ( n = 0; n < nn; n++) {
    fprintf(stderr, "Node %"ISYM"  Memory %8.2f GBytes\n", n, nodemem[n]);
  }
  if ( id == 0 ) fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");

  MPI_Barrier(MPI_COMM_WORLD);

  delete [] grid_size;
  delete [] grid;
  delete [] task_node;
  delete [] task;
  delete [] freemem;
  delete [] freecpu;
  delete [] actual_node;
  delete [] node_input;
  delete [] tasks_per_node;
  delete [] nodemem;

  return SUCCESS;
}

#endif  /* USE_MPI */
#endif  /* ENABLE_TASKMAP */
