#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "nrsrc/nrutil.h"
#include "proto.h"



int comp_func_partminid(void const *a, void const *b)
{
  struct particle_data const *pa, *pb;
  
  pa= a;
  pb= b;

  if(pa->MinID < pb->MinID)
    return -1;

  if(pa->MinID > pb->MinID)
    return +1;
  
  return 0;
}



int comp_func_partcoord(void const *a, void const *b)
{
  struct particle_data const *pa, *pb;
  
  pa= a;
  pb= b;

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


int comp_func(void const *a, void const *b)
{
  struct id_data const *pa, *pb;
  pa= a;
  pb= b;

  if(pa->ID < pb->ID)
    return -1;

  if(pa->ID > pb->ID)
    return +1;

  return 0;
}



int comp_func2(void const *a, void const *b)
{
  struct idmin_data const *pa, *pb;
  pa= a;
  pb= b;

  if(pa->minID < pb->minID)
    return -1;

  if(pa->minID > pb->minID)
    return +1;

  return 0;
}



int comp_func_gr(void const *a, void const *b)
{
  struct gr_data const *pa, *pb;
  pa= a;
  pb= b;

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







