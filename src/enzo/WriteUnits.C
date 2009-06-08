/***********************************************************************
/
/  WRITES UNITS TO AN OUTPUT FILE
/
/  written by: Elizabeth Tasker
/  date:       May, 2005
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "units.h"

int WriteUnits(FILE *fptr)
{

  /* write output to file */

  fprintf(fptr, "MassUnits = %"GSYM"\n", GlobalMassUnits);
  fprintf(fptr, "DensityUnits    = %"GSYM"\n", GlobalDensityUnits);
  fprintf(fptr, "TimeUnits    = %"GSYM"\n", GlobalTimeUnits);
  fprintf(fptr, "LengthUnits  = %"GSYM"\n", GlobalLengthUnits);
     
  return SUCCESS;
}
