/*
   This is just a C wrapper for drand48() called from Fortran code.
*/
 
#include <stdio.h>
#include <stdlib.h>
 
#include "macros_and_parameters.h"

extern "C" FLOAT FORTRAN_NAME(ran1)(int *idim);
 
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
extern "C" FLOAT enzo_ranf_(int *irangen)
#endif
 
#if defined(SP2) || defined(HP)
extern "C" FLOAT enzo_ranf(int *irangen)
#endif
 
{
  int Zero = 0;

  if (*irangen == 0)
    return ((FLOAT) drand48());

  if (*irangen == 1)
    return FORTRAN_NAME(ran1)(&Zero);

  if (*irangen < 0 || *irangen > 1) {
    printf("Value of RandomNumberGenerator = %d is unknown.\n", *irangen);
    exit(EXIT_FAILURE);
  }
}
