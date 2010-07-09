/***********************************************************************
/
/  READS THE RADIATION FIELD DATA
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
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(calc_photo_rates)(
                      int *NFREQ, float *FREQDEL, int *iradshield, float *aye,
		      float *SIGH, float *SIGHE, float *SIGHE2, float *INUTOT,
		      float *PHTH, float *PHTHE, float *PHTHE2,
		          float *EXRY, float *TXRY,
		      float *PHTLAMH, float *PHTLAMHE, float *PHTLAMHE2,
		      float *AVGSIGH, float *AVGSIGHE, float *AVGSIGHE2,
		      float *AVGSIGHP, float *AVGSIGHEP, float *AVGSIGHE2P,
		      float *utim, float *uxyz, float *urho, float *uaye);
 
 
int ReadRadiationData(FILE *fptr)
{
  int i;

 
  /* read in scalar data. */
 
  if (fscanf(fptr, "TimeFieldLastUpdated = %"PSYM,
	     &RadiationData.TimeFieldLastUpdated) != 1) {
    ENZO_FAIL("Error reading TimeFieldLastUpdated.\n");
  }
    
  /* read in field. */

  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) { 
    
    for (i = 0; i < RadiationData.NumberOfFrequencyBins; i++)
      if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM,
		 RadiationData.Spectrum[0]+i, RadiationData.Spectrum[1]+i,
		 RadiationData.Spectrum[2]+i, RadiationData.Spectrum[3]+i)
	  != 4) {
	ENZO_VFAIL("Error reading RadiationData line %"ISYM"\n", i)
      }

  }
 
  /* Compute the units. */
 
  FLOAT a = 1, dadt;
  float DensityUnits = 1, LengthUnits = 1, aUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, 
	       RadiationData.TimeFieldLastUpdated) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(RadiationData.TimeFieldLastUpdated,
					&a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");

    }
    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  float afloat = float(a);
 
  /* Given the field, compute the new radiation-dependant rates. */

  FORTRAN_NAME(calc_photo_rates)(
                &RadiationData.NumberOfFrequencyBins,
		   &RadiationData.FrequencyBinWidth, &RadiationData.RadiationShield, &afloat,
		RadiationData.HICrossSection,
		   RadiationData.HeICrossSection,
		   RadiationData.HeIICrossSection,
		   RadiationData.Spectrum[0],
		&RateData.k24, &RateData.k25, &RateData.k26,
		   &RadiationData.ComptonXrayEnergyDensity,
		   &RadiationData.ComptonXrayTemperature,
		&CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
		&RadiationData.HIAveragePhotoionizationCrossSection,
		   &RadiationData.HeIAveragePhotoionizationCrossSection,
		   &RadiationData.HeIIAveragePhotoionizationCrossSection,
		&RadiationData.HIAveragePhotoHeatingCrossSection,
		   &RadiationData.HeIAveragePhotoHeatingCrossSection,
		   &RadiationData.HeIIAveragePhotoHeatingCrossSection,
		&TimeUnits, &LengthUnits, &DensityUnits, &aUnits);

  return SUCCESS;
}
