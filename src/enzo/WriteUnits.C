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
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "units.h"

int WriteUnits(FILE *fptr)
{

  /* write output to file */

  fprintf(fptr, "MassUnits = %"GOUTSYM"\n", GlobalMassUnits);
  fprintf(fptr, "DensityUnits    = %"GOUTSYM"\n", GlobalDensityUnits);
  fprintf(fptr, "TimeUnits    = %"GOUTSYM"\n", GlobalTimeUnits);
  fprintf(fptr, "LengthUnits  = %"GOUTSYM"\n", GlobalLengthUnits);
     
  return SUCCESS;
}
