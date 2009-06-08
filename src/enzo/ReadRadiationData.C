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


#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
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
    fprintf(stderr, "Error reading TimeFieldLastUpdated.\n");
    return FAIL;
  }
 
  /* read in field. */
 
  for (i = 0; i < RadiationData.NumberOfFrequencyBins; i++)
    if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM,
	       RadiationData.Spectrum[0]+i, RadiationData.Spectrum[1]+i,
	       RadiationData.Spectrum[2]+i, RadiationData.Spectrum[3]+i)
	!= 4) {
      fprintf(stderr, "Error reading RadiationData line %"ISYM"\n", i);
      return FAIL;
    }
 
  /* Compute the units. */
 
  FLOAT a = 1, dadt;
  float DensityUnits = 1, LengthUnits = 1, aUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits,
	       RadiationData.TimeFieldLastUpdated) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(RadiationData.TimeFieldLastUpdated,
					&a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  float afloat = float(a);
 
  /* Given the field, compute the new radiation-dependant rates. */
 
  int RadiationShield = (RadiationFieldType == 10) ? FALSE : TRUE;
  FORTRAN_NAME(calc_photo_rates)(
                &RadiationData.NumberOfFrequencyBins,
		   &RadiationData.FrequencyBinWidth, &RadiationShield, &afloat,
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
