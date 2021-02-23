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

  int Nt; // Temperature
  int Ng; // surface gravity
  int Nz; // metallicity

  float ***q0; // HI photon flux (1/s/cm^2)
  float ***q1; // HeI photon flux (1/s/cm^2)
  float ***q2; // HeII photon flux (1/s/cm^2)

  float ***IR_flux;      // IR  flux (erg/s/cm^2)
  float ***FUV_flux;     // FUV flux (erg/s/cm^2)
  float ***LW_flux;      // LW flux  (erg/s/cm^2)

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

  int Nm;
  int Nz;

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
