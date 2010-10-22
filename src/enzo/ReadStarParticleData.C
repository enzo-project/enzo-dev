/***********************************************************************
/
/  READS THE STAR PARTICLE DATA
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


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
 
 
int ReadStarParticleData(FILE *fptr)
{
 
  if (StarParticleCreation == FALSE)
    return SUCCESS;
 
  /* read in number data. */

  if (!AddParticleAttributes) {
    if (fscanf(fptr, "NumberOfStarParticles = %"ISYM"\n",
	       &NumberOfStarParticles) != 1) {
      //      ENZO_FAIL("Error reading NumberOfStarParticles.\n");

    }
  } else 
    NumberOfStarParticles = 0;

  NumberOfOtherParticles = 0; 

  return SUCCESS;
}
