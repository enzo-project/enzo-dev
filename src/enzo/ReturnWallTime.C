/***********************************************************************
/
/  RETURN THE CURRENT ELAPSED CPU TIME
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

double ReturnWallTime()
{
#ifdef USE_MPI

  return MPI_Wtime();

#else /* USE_MPI */

  return 0.0;

#endif /* USE_MPI */
}
