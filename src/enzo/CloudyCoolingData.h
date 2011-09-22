/***********************************************************************
/
/  Cloudy Cooling Data
/
***********************************************************************/

#define CLOUDY_COOLING_MAX_DIMENSION 5

struct CloudyCoolingDataType
{

  // Use a CMB temperature floor.
  int CMBTemperatureFloor;

  // Flag to control whether or not to include heating from Cloudy.
  int IncludeCloudyHeating;

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  float CloudyElectronFractionFactor;

  // Cooling grid file.
  char *CloudyCoolingGridFile;

  // Rank of Cloudy dataset.
  int CloudyCoolingGridRank;

  // Dimension of Cloudy dataset.
  int *CloudyCoolingGridDimension;

  // Dataset parameter values.
  //  float *CloudyCoolingGridParameters[CLOUDY_COOLING_MAX_DIMENSION];
  float **CloudyCoolingGridParameters;

  // Heating values
  float *CloudyHeating;

  // Cooling values
  float *CloudyCooling;

  // Length of 1D flattened Cloudy data
  int CloudyDataSize;
};
