/***********************************************************************
/
/  WRITES A STAR PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1: Ji-hoon Kim - added additional outputs
/             November, 2009
/  modified2: Michael Kuhlen - moved Ji-hoon's additional output to
/                              [Group_]WriteAllData.C, since this routine
/                              is no longer called when
/                              HierarchyFileInputFormat=0.
/             December 2010
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

    /* mqk 12/04/2010: moved MBH particle data output to
       (Group_)WriteAllData, since this routine is not called when
       HierarchyFileInputFormat=0. */
  }

  return SUCCESS;
}
