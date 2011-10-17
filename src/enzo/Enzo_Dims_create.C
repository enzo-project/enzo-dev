#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

void my_exit(int status);

// External C prototype for MPICH v1.0 routine
extern "C" Eint32 XXMPI_Dims_create(Eint32 nnodes, Eint32 ndims, Eint32 dims[]);


int Enzo_Dims_create(int nnodes, int ndims, int *dims)
{
 
  int nn, mm, i;
  double xn, yn, eps, one_third;

/* Check for cubes */

  if ( ndims == 3 ) {

    one_third = 1.0/3.0;
    eps = 0.1;
    xn = ((double)(nnodes)) + eps;
    yn = POW(xn, one_third);
    nn = (int)(yn);
    mm = nn * nn * nn;

    if ( mm == nnodes ) {
      for ( i = 0; i < ndims; i++ ) {
        dims[i] = nn;
      }
      return SUCCESS;
    } 

  }

/* Not 3D and cubic */

  MPI_Arg mcpu, rank;
  MPI_Arg mpi_layout[] = {0, 0, 0};;

  mcpu = nnodes;
  rank = ndims;

  for ( i = 0; i < 3; i++ ) {
    mpi_layout[i] = 0;
  }

  XXMPI_Dims_create(mcpu, rank, mpi_layout);

  for ( i = 0; i < ndims; i++ ) {
    dims[i] = mpi_layout[i];
  }

  return SUCCESS;
}
