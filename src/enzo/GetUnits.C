/***********************************************************************
/
/  COMPUTE AND RETURN THE UNITS
/
/  written by: Elizabeth Tasker
/  date:       May, 2005
/  modified1:
/
/  PURPOSE:  Returns the units in CGS:
/
/  
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "units.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "TopGridData.h"

/* function prototypes */
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, FLOAT Time);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time)

{
  /* If using cosmology, get cosmology units */
  if (ComovingCoordinates) {
    //fprintf(stderr, "Using CosmologyCoordinates.\n");
    if (CosmologyGetUnits(DensityUnits, LengthUnits, TemperatureUnits,
                          TimeUnits, VelocityUnits, Time) == FAIL) {
      ENZO_FAIL("Error in CosmologyGetUnits.\n");
    }
  }
  else {
      /* Determine the units. */
      *DensityUnits = GlobalDensityUnits;
      *MassUnits = GlobalMassUnits;
      *LengthUnits      = GlobalLengthUnits;
      *TemperatureUnits = mh*pow(GlobalLengthUnits/GlobalTimeUnits,2)/kboltz;   //K
      *TimeUnits        = GlobalTimeUnits;
      *VelocityUnits    = GlobalLengthUnits/GlobalTimeUnits; //cms-1
  
    }
  return SUCCESS;
}

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time)

{
  /* If using cosmology, get cosmology units */
  if (ComovingCoordinates) {
    //    fprintf(stderr, "Using CosmologyCoordinates.\n");
    if (CosmologyGetUnits(DensityUnits, LengthUnits, TemperatureUnits,

                          TimeUnits, VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
  }
  else {
      /* Determine the units. */
      *DensityUnits = GlobalDensityUnits;
      *LengthUnits      = GlobalLengthUnits;
      *TemperatureUnits = mh*pow(GlobalLengthUnits/GlobalTimeUnits,2)/kboltz;   //K
      *TimeUnits        = GlobalTimeUnits;
      *VelocityUnits    = GlobalLengthUnits/GlobalTimeUnits; //cms-1
    }
  return SUCCESS;
}
