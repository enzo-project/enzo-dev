/***********************************************************************
/
/  WRITES A STAR PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1: Ji-hoon Kim - added additional outputs
/             November, 2009
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
#include "TopGridData.h"
#include "StarParticleData.h"
 
 
int WriteStarParticleData(FILE *fptr, TopGridData &MetaData)
{
 
  if (StarParticleCreation == FALSE)
    return SUCCESS;

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Write out particle number data. */

    fprintf(fptr, "\n");    
    fprintf(fptr, "NumberOfStarParticles      = %"ISYM"\n", NumberOfStarParticles);
    fprintf(fptr, "NumberOfOtherParticles     = %"ISYM"\n", NumberOfOtherParticles); 

  }

  return SUCCESS;
}
