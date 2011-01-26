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
 
  if (StarParticleCreation == FALSE && MBHParticleIO == FALSE)
    return SUCCESS;

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Write out particle number data. */

    fprintf(fptr, "\n");    
    fprintf(fptr, "NumberOfStarParticles      = %"ISYM"\n", NumberOfStarParticles);
    fprintf(fptr, "NumberOfOtherParticles     = %"ISYM"\n", NumberOfOtherParticles); 

    /* Write out MBH particle data (mass, angular momentum) */ 

    if(MBHParticleIO == TRUE) {

      FILE *fptr2;
      if ((fptr2 = fopen(MBHParticleIOFilename, "a")) == NULL) {
	ENZO_VFAIL("Error opening file %s\n", MBHParticleIOFilename)

      }

      // printing order: time, regular star count, MBH id, MBH mass, MBH angular momentum
      for (int i = 0; i < G_TotalNumberOfStars; i++) { 
	fprintf(fptr2, " %"FSYM"  %"ISYM"  %"ISYM"  %lf  %"FSYM"  %"FSYM"  %"FSYM"  %lf\n", 
		MetaData.Time, NumberOfStarParticles, (int)(MBHParticleIOTemp[i][0]), 
		MBHParticleIOTemp[i][1], (float)(MBHParticleIOTemp[i][2]), 
		(float)(MBHParticleIOTemp[i][3]), (float)(MBHParticleIOTemp[i][4]),
		MBHParticleIOTemp[i][5]);
      }

      fclose(fptr2);
    }

  }

  return SUCCESS;
}
