/***********************************************************************
/
/  Cosmic Ray Acceleration Efficiency Data		
/
***********************************************************************/

struct CosmicRayDataType
{
  // Number of dimension to interpolate over
  int CosmicRayGridRank;

  // Maximum number of values for Pre-CR population
  int CRMaxNumberPrePopValues;
  // Maximum number of values for Mach 
  int CRMaxNumberMachValues;

  // Cooling grid run file
  char *CosmicRayGridRunFile;

  // Values for first parameter space (Pre-Shock CR Percentage)
  float *CRPrePopValues;
  // Values for second parameter space (Mach Number)
  float *CRMachValues;

  //Min/max of the paramters
  float CRMinPrePop;
  float CRMaxPrePop;
  float CRMinMach;
  float CRMaxMach;

  // Number of values for first parameter space
  int CRNumberPrePopValues;
  // Number of values for second parameter space
  int CRNumberMachValues;

  // Efficiency Values
  float **CREfficiency;
};
