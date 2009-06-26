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
  

  /*  metallicity fields */
  int UseMetallicityField;
  float MetallicityField_Fraction;

  int MultiMetals;
  float MultiMetalsField1_Fraction;
  float MultiMetalsField2_Fraction;

  /* Cooling Test parameters */
  float MinimumHNumberDensity;
  float MaximumHNumberDensity;
  float MinimumMetallicity;
  float MaximumMetallicity;
  float MinimumTemperature;
  float MaximumTemperature;
  int ResetEnergies;

};
