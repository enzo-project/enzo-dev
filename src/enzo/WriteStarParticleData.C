/***********************************************************************
/
/  WRITES A STAR PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>


#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
 
 
int WriteStarParticleData(FILE *fptr)
{
 
  if (StarParticleCreation == FALSE)
    return SUCCESS;
 
  /* Write out number data. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfStarParticles      = %"ISYM"\n", NumberOfStarParticles);
 
  return SUCCESS;
}
