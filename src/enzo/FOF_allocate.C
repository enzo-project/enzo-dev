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


void deallocate_all_memory(FOFData &D)
{

  delete [] D.Nslab;
  delete [] D.Nshadow;
  delete [] D.Noffset;
  delete [] D.NtoLeft;
  delete [] D.NtoRight;
  delete [] D.ContribID;
  delete [] D.ContribHead;
  delete [] D.GridNext;
  delete [] D.Head;
  delete [] D.Tail;
  delete [] D.Next;
  delete [] D.Len;
  delete [] D.GroupDat;
  delete [] D.GroupDatAll;
  delete [] D.P;

  return;

}
