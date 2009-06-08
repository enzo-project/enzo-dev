/***********************************************************************
/
/  MULTI-SPECIES RATE DATA (see solve_rate.src)
/
***********************************************************************/

struct RateDataType
{
  int NumberOfTemperatureBins;   
  float TemperatureStart;        // range of temperature in K
  float TemperatureEnd;

  /* 6 species rates */

  float *k1;
  float *k2;
  float *k3;
  float *k4;
  float *k5;
  float *k6;

  /* 9 species rates (including H2) */

  float *k7;
  float *k8;
  float *k9;
  float *k10;
  float *k11;
  float *k12;
  float *k13;
  float *k14;
  float *k15;
  float *k16;
  float *k17;
  float *k18;
  float *k19;
  float *k20;  /* currently not used */
  float *k21;  /* currently not used */
  float *k22;  /* 3-body H2 formation */

  float *k13dd;  /* density dependent version of k13 (collisional H2
                    dissociation); actually 7 functions instead of 1. */

  /* Radiative rates for 6-species (for external field). */

  float k24;
  float k25;
  float k26;

  /* Radiative rates for 9-species (for external field). */

  float k27;
  float k28;
  float k29;
  float k30;
  float k31;

  /* 12 species rates (with Deuterium). */

  float *k50;
  float *k51;
  float *k52;
  float *k53;
  float *k54;
  float *k55;
  float *k56;
};
