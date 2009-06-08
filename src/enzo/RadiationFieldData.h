/***********************************************************************
/
/  RADIATION FIELD DATA
/
***********************************************************************/

struct RadiationFieldDataType
{

  FLOAT TimeFieldLastUpdated;
  int  RadiationShield;

  /* Frequency dependant information (constant). */

  int NumberOfFrequencyBins;
  float FrequencyBinWidth;    // in log10(eV)

  float *HICrossSection;
  float *HeICrossSection;
  float *HeIICrossSection;

  float *StellarSpectrum;
  float *QuasarSpectrum;

  /* Frequency dependant information (current spectra)
      0 - Total, 1 - Gas only, 2 - stars only, 3 - Quasar only. */

  float *Spectrum[4];
  float *Emissivity[4];

  /* Temperature dependant information. */

  int NumberOfTemperatureBins;
  float TemperatureBinWidth;
  float TemperatureBinMinimum;

  float *FreeFreeDensity[MAX_DEPTH_OF_HIERARCHY];
  float *HIIFreeBoundDensity[MAX_DEPTH_OF_HIERARCHY];
  float *HeIIFreeBoundDensity[MAX_DEPTH_OF_HIERARCHY];
  float *HeIIIFreeBoundDensity[MAX_DEPTH_OF_HIERARCHY];

  /* global averages */

  float HIMeanDensity[MAX_DEPTH_OF_HIERARCHY];
  float HeIMeanDensity[MAX_DEPTH_OF_HIERARCHY];
  float HeIIMeanDensity[MAX_DEPTH_OF_HIERARCHY];

  float HIAveragePhotoionizationCrossSection;
  float HeIAveragePhotoionizationCrossSection;
  float HeIIAveragePhotoionizationCrossSection;
  float HIAveragePhotoHeatingCrossSection;
  float HeIAveragePhotoHeatingCrossSection;
  float HeIIAveragePhotoHeatingCrossSection;

  float ComptonXrayTemperature;
  float ComptonXrayEnergyDensity;

  float IntegratedStarFormation;

};
