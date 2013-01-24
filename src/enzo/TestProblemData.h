/***********************************************************************
/
/  Test problem data structure
/  Look in SetDefaultGlobalValues.C for default settings!
/
***********************************************************************/
struct TestProblemDataType
{

  float HydrogenFractionByMass;
  float DeuteriumToHydrogenRatio;

  /* multispecies */
  int MultiSpecies;

  // fractions in the interior (where the blast wave energy is injected)
  float HI_Fraction_Inner;
  float HII_Fraction_Inner;
  float HeI_Fraction_Inner;
  float HeII_Fraction_Inner;
  float HeIII_Fraction_Inner;
  float HM_Fraction_Inner;
  float H2I_Fraction_Inner;
  float H2II_Fraction_Inner;
  float DI_Fraction_Inner;
  float DII_Fraction_Inner;
  float HDI_Fraction_Inner;

  // fractions in the surrouding medium
  float HI_Fraction;
  float HII_Fraction;
  float HeI_Fraction;
  float HeII_Fraction;
  float HeIII_Fraction;
  float HM_Fraction;
  float H2I_Fraction;
  float H2II_Fraction;
  float DI_Fraction;
  float DII_Fraction;
  float HDI_Fraction;
  
  int UseMassInjection;
  float InitialHydrogenMass;
  float InitialHeliumMass;
  float InitialDeuteriumMass;

  // H mass fraction for inner region in blast wave problem
  float InnerHydrogenFractionByMass;
  // D/H for inner region in blast wave problem
  float InnerDeuteriumToHydrogenRatio;

  /*  metallicity fields */
  int UseMetallicityField;
  float MetallicityField_Fraction;
  float MetallicitySNIaField_Fraction;
  float MetallicityNormalization;

  float InitialMetalMass;

  int MultiMetals;
  float MultiMetalsField1_Fraction;
  float MultiMetalsField2_Fraction;

  /* simon glover chemistry+cooling */
  int GloverChemistryModel;
  
  float StartOfHighTempRegime;

  // fractions in the interior (where the blast wave energy is injected)
  float CI_Fraction_Inner;
  float CII_Fraction_Inner;
  float OI_Fraction_Inner;
  float OII_Fraction_Inner;
  float SiI_Fraction_Inner;
  float SiII_Fraction_Inner;
  float SiIII_Fraction_Inner;
  float CHI_Fraction_Inner;
  float CH2I_Fraction_Inner;
  float CH3II_Fraction_Inner;
  float C2I_Fraction_Inner;
  float COI_Fraction_Inner;
  float HCOII_Fraction_Inner;
  float OHI_Fraction_Inner;
  float H2OI_Fraction_Inner;
  float O2I_Fraction_Inner;

  // fractions in the surrouding medium
  // these are also the fractions used in the uniform grid setup 
  // and constant density problem.
  float CI_Fraction;
  float CII_Fraction;
  float OI_Fraction;
  float OII_Fraction;
  float SiI_Fraction;
  float SiII_Fraction;
  float SiIII_Fraction;
  float CHI_Fraction;
  float CH2I_Fraction;
  float CH3II_Fraction;
  float C2I_Fraction;
  float COI_Fraction;
  float HCOII_Fraction;
  float OHI_Fraction;
  float H2OI_Fraction;
  float O2I_Fraction;

  /* Cooling Test parameters */
  float MinimumHNumberDensity;
  float MaximumHNumberDensity;
  float MinimumMetallicity;
  float MaximumMetallicity;
  float MinimumTemperature;
  float MaximumTemperature;
  int ResetEnergies;

  /* Shock Fields */
  int ShockMethod;
  int StorePreShockFields;
  float ShockTemperatureFloor;

  /* constant for analytical solution to free-fall collapse */
  float OneZoneFreefallConstant;
  /* fraction of free-fall time for timestep */
  float OneZoneFreefallTimestepFraction;
  float OneZoneFreefallUseEffectiveGamma;
};
