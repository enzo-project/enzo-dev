/***********************************************************************
/
/  INITIALIZE THE RADIATION FIELD
/
/  written by: Greg Bryan
/  date:       October, 1999
/  modified1:
/
/  PURPOSE:
/    For runs with a dynamic internal radiation field, initialize the data
/    (This applies only for RadiationFieldType = 10 or 11)
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
 
extern "C" void FORTRAN_NAME(chtable)(int *NFREQ, float *FREQDEL,
				      float *SIGH, float *SIGHE, float *SIGHE2,
				      float *SPSTAR1, float *SPSTAR2);
 
int InitializeRadiationFieldData(FLOAT Time)
{
 
  int i, j, level;
 
  /* Set radiation data parameters. */
  /* These should really be read in. */
 
  int nfbins = 400, ntbins = 600;
 
  RadiationData.TimeFieldLastUpdated = Time;
 
  RadiationData.NumberOfFrequencyBins = nfbins;
  RadiationData.FrequencyBinWidth = 0.02;    // in log10(eV)
 
  RadiationData.NumberOfTemperatureBins = ntbins;
  RadiationData.TemperatureBinWidth = 0.02;  // in log10(K)
  RadiationData.TemperatureBinMinimum = 3.0;  // in log10(K)
 
  /* Allocate space. */
 
  RadiationData.HICrossSection = new float[nfbins];
  RadiationData.HeICrossSection = new float[nfbins];
  RadiationData.HeIICrossSection = new float[nfbins];
 
  RadiationData.StellarSpectrum = new float[nfbins];
  RadiationData.QuasarSpectrum = new float[nfbins];
 
  for (i = 0; i < 4; i++) {
    RadiationData.Spectrum[i] = new float[nfbins];
    RadiationData.Emissivity[i] = new float[nfbins];
  }
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    RadiationData.FreeFreeDensity[level] = new float[ntbins];
    RadiationData.HIIFreeBoundDensity[level] = new float[ntbins];
    RadiationData.HeIIFreeBoundDensity[level] = new float[ntbins];
    RadiationData.HeIIIFreeBoundDensity[level] = new float[ntbins];
  }
 
  /* Zero fields. */
 
  for (i = 0; i < 4; i++)
    for (j = 0; j < nfbins; j++) {
      RadiationData.Spectrum[i][j] = 0;
      RadiationData.Emissivity[i][j] = 0;
    }
 
  RadiationData.HIAveragePhotoionizationCrossSection = 0; 
  RadiationData.HeIAveragePhotoionizationCrossSection = 0;
  RadiationData.HeIIAveragePhotoionizationCrossSection = 0;
  RadiationData.HIAveragePhotoHeatingCrossSection = 0;
  RadiationData.HeIAveragePhotoHeatingCrossSection = 0;
  RadiationData.HeIIAveragePhotoHeatingCrossSection = 0;
 
  RadiationData.IntegratedStarFormation = 0;
 
  RadiationData.ComptonXrayTemperature = 0;
  RadiationData.ComptonXrayEnergyDensity = 0;
 
  /* Call fortran routine to generate tables. */
 
   FORTRAN_NAME(chtable)(&RadiationData.NumberOfFrequencyBins,
 			&RadiationData.FrequencyBinWidth,
 			RadiationData.HICrossSection,
 			RadiationData.HeICrossSection,
 			RadiationData.HeIICrossSection,
 			RadiationData.StellarSpectrum,
 			RadiationData.QuasarSpectrum);
 
  return SUCCESS;
}
