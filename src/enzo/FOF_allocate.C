#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_proto.h"

void allocate_memory(FOFData &AllVars)
{
  int n;

  n = AllVars.Nslab[MyProcessorNumber] + AllVars.Nshadow[MyProcessorNumber];

  if (n > 0) {
    AllVars.P = new FOF_particle_data[n];
    if (AllVars.P == NULL)
      ENZO_FAIL("failed to allocate memory. (A)");
    AllVars.P -= 1;
  }
}









