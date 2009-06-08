/*
   This is just a C wrapper for srand48() called from Fortran code.
*/
 
#include <stdio.h>
#include <stdlib.h>
 
#include "macros_and_parameters.h"
 
/*
  amd    trailing underscore_            nec    trailing underscore_
  dec    trailing underscore_            sgi    trailing underscore_
  intel  trailing underscore_            sun    trailing underscore_
  linux  trailing underscore_
 
  hp     NO trailing underscore          cray   NO trailing underscore
  ibm    NO trailing underscore
 
  Cray uses UPPERCASE names.
*/
 
#if defined(IRIS4) || defined(SUN)   || defined(COMPAQ) || \
    defined(IA64) || defined(LINUX) || defined(NEC) || defined(CRAYX1) || defined(XT3)
extern "C" void enzo_seed_(int *seed)
#endif
 
#if defined(SP2) || defined(HP)
extern "C" void enzo_seed(int *seed)
#endif
 
{
  long seed_value;
  seed_value = (long) (*seed);
  printf("Seed value %ld\n", seed_value);
  srand48(seed_value);
}
