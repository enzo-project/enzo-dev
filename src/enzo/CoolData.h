/***********************************************************************
/
/  MULTI-SPECIES COOLING DATA (see multi_cool.src)
/
***********************************************************************/

struct CoolDataType
{
  int NumberOfTemperatureBins;   
  int ih2co;                     // flag for H2 cooling (0-off/1-on)
  int ipiht;                     // flag for photoionization cooling
  float TemperatureStart;        // range of temperature in K
  float TemperatureEnd;

  /* Radiation parameters (should be put in own spot). */

  float alpha0;
  float f3;
  float f0to3;
  float RadiationRedshiftOn;
  float RadiationRedshiftOff;
  float RadiationRedshiftFullOn;
  float RadiationRedshiftDropOff;
  float HydrogenFractionByMass;
  float DeuteriumToHydrogenRatio;
  float SolarMetalFractionByMass;

  /* Equilibrium rates */

  float *EquilibriumRate;

  /* 6 species rates */

  float *ceHI;                   // collisional excitation rates
  float *ceHeI;
  float *ceHeII;
  float *ciHI;                   // collisional ionization
  float *ciHeI;
  float *ciHeIS;
  float *ciHeII;
  float *reHII;                  // recombination
  float *reHeII1;
  float *reHeII2;
  float *reHeIII;
  float *brem;                   // free-free (Bremsstrahlung)
  float comp;                    // Compton cooling
  float comp_xray;               // X-ray compton heating coefficient
  float temp_xray;               // X-ray compton heating temperature (K)
  float gammah;                  // Photoelectric heating

  /* radiative rates (external field). */

  float piHI;                    // photo-ionization cooling
  float piHeI;                   //    (no temperature dependance)
  float piHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi.src. */

  float *hyd01k;
  float *h2k01;
  float *vibh;
  float *roth;
  float *rotl;

  float *GP99LowDensityLimit;
  float *GP99HighDensityLimit;

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  float *GAHI;
  float *GAH2;
  float *GAHe;
  float *GAHp;
  float *GAel;

  /* 12 species rates (including HD) */

  float *HDlte;
  float *HDlow;
  float *HDcool;

  /* CIE cooling */
  float *cieco;

  /* Metal line cooling (tabulated by temperature and electron
     fraction) */

  int NumberOfElectronFracBins;
  float ElectronFracStart;
  float ElectronFracEnd;
  float *metals;

  /* Energy transfer to grains. */
  float *gas_grain;

  /* For analysis, ratios of metal fine structure line emission is
     desired sometimes. */

  float *metal_ratios;
  int MR_NumberOfTemperatureBins;
  int MR_NumberOfElectronFracBins;
  float MR_TemperatureStart;
  float MR_TemperatureEnd;
  float MR_ElectronFracStart;
  float MR_ElectronFracEnd;

};
