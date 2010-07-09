/***********************************************************************
/
/  INITIALIZE THE SPECTRUM TABLE FOR RADIATIVE TRANSFER
/
/  written by: Ji-hoon Kim
/  date:       February, 2010
/  modified1:
/
/  PURPOSE: For RadiativeTransferTraceSpectrum = TRUE, 
/           read in and initialize the spectrum table.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#ifdef TRANSFER
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int ReadRadiativeTransferSpectrumTable(float TemperatureUnits, float LengthUnits, 
				       float aUnits, float DensityUnits, float TimeUnits);
 
int InitializeRadiativeTransferSpectrumTable(FLOAT Time)
{
 
  FLOAT a = 1, dadt;
 
  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }

  /* Read spectrum table */

  if (RadiativeTransferTraceSpectrum == TRUE) {
    if (ReadRadiativeTransferSpectrumTable(TemperatureUnits, LengthUnits, aUnits, 
					   DensityUnits, TimeUnits) == FAIL) {
      ENZO_FAIL("Error in ReadRadiativeTransferSpectrumTable.\n");

    }
  }

  return SUCCESS;
}

#endif
