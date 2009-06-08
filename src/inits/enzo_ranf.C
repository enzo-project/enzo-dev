/*
   This is just a C wrapper for drand48() called from Fortran code.
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
extern "C" FLOAT enzo_ranf_(void)
#endif
 
#if defined(SP2) || defined(HP)
extern "C" FLOAT enzo_ranf(void)
#endif
 
{
  return ((FLOAT) drand48());
}
