/***********************************************************************
/
/  WRITES THE RADIATION FIELD DATA
/
/  written by: Greg Bryan
/  date:       October, 1999
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
 
 
int WriteRadiationData(FILE *fptr)
{
  int i;
 
  /* write scalar data. */
 
  fprintf(fptr, "TimeFieldLastUpdated = %"GOUTSYM"\n",
	  RadiationData.TimeFieldLastUpdated);
 
  /* write field. */

  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) {

    for (i = 0; i < RadiationData.NumberOfFrequencyBins; i++)
      fprintf(fptr, "%"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	      RadiationData.Spectrum[0][i], RadiationData.Spectrum[1][i],
	      RadiationData.Spectrum[2][i], RadiationData.Spectrum[3][i]);
    
  }
 
  return SUCCESS;
}
