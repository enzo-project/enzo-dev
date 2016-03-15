/*-------------------------------------------------------------------------
/
/ INDIVIDUAL STAR DATA
/
/ Author: Andrew Emerick
/ Date  : March, 2016
/ modified1:
/
/ PURPOSE:
/  Definition of struct to store radiation data look up table for
/  computing ionizing fluxes of individual stars
/
---------------------------------------------------------------------------*/

struct IndividualStarRadDataType
{

  int NumberOfTemperatureBins;
  int NumberOfSGBins;
  int NumberOfMetallicityBins;

//  float **q0[INDIVIDUAL_STAR_TEMPERATURE_BINS];
//  float **q1[INDIVIDUAL_STAR_TEMPERATURE_BINS];
  float ***q0;
  float ***q1;

  // bin values
  float *T;
  float *Z;
  float *g;
};

