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

#include "ErrorExceptions.h"
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
    ret += sscanf(line, "MassUnits = %lf", &GlobalMassUnits);
    ret += sscanf(line, "DensityUnits = %"FSYM, &GlobalDensityUnits);
    ret += sscanf(line, "LengthUnits = %"FSYM, &GlobalLengthUnits);
    ret += sscanf(line, "TimeUnits = %"FSYM, &GlobalTimeUnits);
  }

  /* If both mass and density units specified, use only mass and print warning */

  if (GlobalDensityUnits != 1.0 && GlobalMassUnits != 1.0) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Warning! Density and Mass units are both defined. "
	      "Using only MassUnits.\n");
    GlobalDensityUnits = 1.0;
  }

  /* if only one of density/mass units specifed, calculate the other one */
  if (GlobalMassUnits == 1.0)
    GlobalMassUnits = (double) (double(GlobalDensityUnits) * pow(GlobalLengthUnits,3));
  else if (GlobalDensityUnits == 1.0)
    GlobalDensityUnits = (float) (double(GlobalMassUnits) / pow(GlobalLengthUnits,3));

  /* When Gravity is used one does not have the freedom to set Length, Density and time units */
  /* We blindly assume here  that if you specified the DensityUnits you want to set the time units
     accordingly.  Tom Abel 2009  
  I doubt this will cause a problem ... but when you read this it probably did ... */
  if (GlobalTimeUnits == 1 && GlobalDensityUnits != 1) {
    if(SelfGravity) {
      GlobalTimeUnits =  1/sqrt(6.67428e-8*GlobalDensityUnits);
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "****** ReadUnits: Set Time Units based on Density Units "
		"u_t = 1./sqrt(G u_rho).\n");
    }
    else {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "****** Warning! TimeUnits = 1.0 at the moment! \n");      
	fprintf(stderr, "****** You have to choose your own TimeUnits when SelfGravity = 0.\n");      
      }
    }
  }

  if (debug)
    fprintf(stderr,"****** ReadUnits:  %e %e %e %e *******\n",
	    GlobalMassUnits, GlobalDensityUnits, GlobalLengthUnits,
	    GlobalTimeUnits);

  return SUCCESS;
}
