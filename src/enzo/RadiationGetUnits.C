/***********************************************************************
/
/  COMPUTE AND RETURN THE RADIATION UNITS
/
/  written by: Daniel R. Reynolds
/  date:       June 2008
/  modified1:
/
/  PURPOSE:  Returns the radiation units in CGS:
/
/  
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "units.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"

/* function prototypes */


int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int RadiationGetUnits(float *RadiationUnits, FLOAT Time)

{

  double MassUnits=1.0;
  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, 
    TimeUnits=1.0, VelocityUnits=1.0;

  // Get Enzo units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  
  *RadiationUnits = DensityUnits*VelocityUnits*VelocityUnits;

  return SUCCESS;
}
