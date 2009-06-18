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
#include "FOF_nrutil.h"
#include "FOF_proto.h"


Eint32 comp_func_partminid(void const *a, void const *b)
{
  FOF_particle_data *pa, *pb;
  
  pa = (FOF_particle_data*) a;
  pb = (FOF_particle_data*) b;

  if(pa->MinID < pb->MinID)
    return -1;

  if(pa->MinID > pb->MinID)
    return +1;
  
  return 0;
}



Eint32 comp_func_partcoord(void const *a, void const *b)
{
  FOF_particle_data *pa, *pb;
  
  pa = (FOF_particle_data*) a;
  pb = (FOF_particle_data*) b;

  if(pa->Pos[0]<pb->Pos[0])
    return -1;

  if(pa->Pos[0]>pb->Pos[0])
    return +1;

  if(pa->Pos[1]<pb->Pos[1])
    return -1;

  if(pa->Pos[1]>pb->Pos[1])
    return +1;

  if(pa->Pos[2]<pb->Pos[2])
    return -1;

  if(pa->Pos[2]>pb->Pos[2])
    return +1;

  return 0;
}


Eint32 comp_func(void const *a, void const *b)
{
  id_data *pa, *pb;
  pa = (id_data*) a;
  pb = (id_data*) b;

  if(pa->ID < pb->ID)
    return -1;

  if(pa->ID > pb->ID)
    return +1;

  return 0;
}



Eint32 comp_func2(void const *a, void const *b)
{
  idmin_data *pa, *pb;
  pa = (idmin_data*) a;
  pb = (idmin_data*) b;

  if(pa->minID < pb->minID)
    return -1;

  if(pa->minID > pb->minID)
    return +1;

  return 0;
}



Eint32 comp_func_gr(void const *a, void const *b)
{
  gr_data *pa, *pb;
  pa = (gr_data*) a;
  pb = (gr_data*) b;

  if(pa->Len < pb->Len)
    return -1;

  if(pa->Len > pb->Len)
    return +1;

  if(pa->Tag < pb->Tag)
    return -1;

  if(pa->Tag > pb->Tag)
    return +1;

  return 0;
}

Eint32 compare_slab(const void *a, const void *b)
{
  struct FOF_particle_data *ia = (struct FOF_particle_data*) a;
  struct FOF_particle_data *ib = (struct FOF_particle_data*) b;
  if (ia->slab - ib->slab < 0)
    return -1;
  else if (ia->slab - ib->slab > 0)
    return 1;
  return 0;
}
