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

  n = AllVars.Nslab[MyProcessorNumber] + AllVars.Nshadow[MyProcessorNumber] + 1;

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
  delete [] D.GroupDat;
  delete [] D.GroupDatAll;

  free_ivector(D.ContribID, 0, D.Ncontrib-1);
  free_ivector(D.ContribHead, 0, D.Ncontrib-1);
  if (D.Nlocal > 0) {
    free_ivector(D.Head, 1, D.Nlocal);
    free_ivector(D.Tail, 1, D.Nlocal);
    free_ivector(D.Next, 1, D.Nlocal);
    free_ivector(D.Len, 1, D.Nlocal);
    free_ivector(D.GridNext, 1, D.Nlocal);
  }

  free_i3tensor(D.GridFirst, 0, D.Grid-1, 0, D.Grid-1, 0, D.Grid-1);
  free_i3tensor(D.GridLast,  0, D.Grid-1, 0, D.Grid-1, 0, D.Grid-1);
  free_i3tensor(D.GridFlag,  0, D.Grid-1, 0, D.Grid-1, 0, D.Grid-1);

  /* Now we do this after copying back to grids in FOF_Finalize */

//  D.P++;
//  delete [] D.P;

  return;

}
