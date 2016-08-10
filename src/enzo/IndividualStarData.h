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

  float ***Fuv; // FUV heating rate

  // bin values
  float *T;
  float *Z;
  float *g;

  // solar metallicity as defined by OSTAR2002
  // using Grevesse & Savaul 1998 0.017
  float Zsolar;
};

struct IndividualStarPropertiesDataType
{

  int NumberOfMassBins;
  int NumberOfMetallicityBins;

  // interpolation values
  float **Teff;
  float **R;
  float **L;
  float **lifetime;
  float **agb_start;

  // bin values
  float *M; // mass (solar)
  float *Z; // metal fraction (NOT SOLAR UNITS)

  // solar metallicity as defined in PARSEC
  // stellar evolution code using
  // Caffau et. al. 2009/2011  0.01524
  float Zsolar;
};
