/***********************************************************************
/
/  READS UNITS FROM INPUT FILE
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
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"

int ReadUnits(FILE *fptr)
{

  char line[MAX_LINE_LENGTH];

  /* Set defaults. */

  GlobalMassUnits        = 1.0;
  GlobalDensityUnits     = 1.0;
  GlobalLengthUnits      = 1.0;
  GlobalTimeUnits        = 1.0;
  
  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;

    /* read parameters */
    
    ret += sscanf(line, "MassUnits = %"FSYM, &GlobalMassUnits);
    ret += sscanf(line, "DensityUnits = %"FSYM, &GlobalDensityUnits);
    ret += sscanf(line, "LengthUnits = %"FSYM, &GlobalLengthUnits);
    ret += sscanf(line, "TimeUnits = %"FSYM, &GlobalTimeUnits);
  }
      

  /* If both mass and density units specified, use only mass and print warning */

  if (GlobalDensityUnits != 1.0 && GlobalMassUnits != 1.0){
    fprintf(stderr, "Warning! Density and mass units are both defined. Using only mass units\n");
    GlobalDensityUnits = 1.0;
  }


  /* if only one of density/mass units specifed, calculate the other one */

  if (GlobalMassUnits == 1.0)
    GlobalMassUnits = (float) (double(GlobalDensityUnits) * pow(GlobalLengthUnits,3));

  else if (GlobalDensityUnits == 1.0)
    GlobalDensityUnits = (float) (double(GlobalMassUnits) / pow(GlobalLengthUnits,3));


  fprintf(stderr,"****** GetUnits:  %e %e %e %e *******\n",GlobalMassUnits,GlobalDensityUnits,GlobalLengthUnits, GlobalTimeUnits);

  return SUCCESS;
}
