/***********************************************************************
/
/  INITIALIZE THE MULTI-SPECIES RATES
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  Dan Reynolds, July 2010; added case-B recombination rates
/  modified2:  Britton Smith, October 2010; moved reading/writing of 
/              parameters to Read/WriteParameterFile.
/
/  PURPOSE:
/    For multi-species runs (with cooling), initialize both the
/      CoolData and RateData rate tables.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int InitializeCloudyCooling(FLOAT Time);
int InitializeLymanWernerTable();
int InitializeHM12Photorates();
int ReadMetalCoolingRates(float TemperatureUnits, float LengthUnits, 
			  float aUnits, float DensityUnits, float TimeUnits, 
			  float aye);
extern "C" void FORTRAN_NAME(calc_rates)(
     int *nratec, float *aye, float *temstart, float *temend, float *alpha0,
     float *f3, int *iradtype, int *casebrates, int *threebody,
     float *utem, float *uxyz, float *uaye, float *urho, float *utim,
     float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
     float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
     float *reHeII2a, float *reHeIIIa, float *brema, float *compa, 
     float *gammahacgs, float *gammaha,
     float *piHI, float *piHeI, float *piHeII,
     float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
     float *gpldl, float *gphdl, float *hdlte, float *hdlow, float *hdcool, float *cieco,
     float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela, float *gasgr, 
     float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, float *k6a,
        float *k7a, float *k8a, float *k9a, float *k10a,
     float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a,
        float *k15a, float *k16a, float *k17a, float *k18a,
     float *k19a, float *k20a, float *k21a, float *k22, float *k23,
     float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
        float *k30, float *k31,
     float *k50, float *k51, float *k52, float *k53, float *k54, float *k55,
        float *k56, int *ndratec, float *dtemstart, float *dtemend, float *h2dusta, 
     float *ncrca, float *ncrd1a, float *ncrd2a, int *ioutput);


// character strings
EXTERN char outfilename[];

 
int InitializeRateData(FLOAT Time)
{
 
  /* Declarations. */
 
  FLOAT a = 1, dadt;

  if (H2FormationOnDust && !RadiativeCooling) {
    ENZO_FAIL("For H2FormationOnDust = 1, must have RadiativeCooling = 1 and MultiSpecies > 0.");
  }
  
  if (debug) printf("InitializeRateData: NumberOfTemperatureBins = %"ISYM"\n",
		    CoolData.NumberOfTemperatureBins);

  if (debug) printf("InitializeRateData: RadiationFieldType = %"ISYM"\n",
		    RadiationFieldType);


  /* Allocate CoolData space for rates. */
 
  CoolData.ceHI    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ceHeI   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ceHeII  = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ciHI    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ciHeI   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ciHeIS  = new float[CoolData.NumberOfTemperatureBins];
  CoolData.ciHeII  = new float[CoolData.NumberOfTemperatureBins];
  CoolData.reHII   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.reHeII1 = new float[CoolData.NumberOfTemperatureBins];
  CoolData.reHeII2 = new float[CoolData.NumberOfTemperatureBins];
  CoolData.reHeIII = new float[CoolData.NumberOfTemperatureBins];
  CoolData.brem    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.hyd01k  = new float[CoolData.NumberOfTemperatureBins];
  CoolData.h2k01   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.vibh    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.roth    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.rotl    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GP99LowDensityLimit  = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GP99HighDensityLimit = new float[CoolData.NumberOfTemperatureBins];
  CoolData.HDlte   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.HDlow   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.HDcool  = new float[CoolData.NumberOfTemperatureBins*5];
  CoolData.cieco   = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAHI    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAH2    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAHe    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAHp    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAel    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.gas_grain = new float[CoolData.NumberOfTemperatureBins]; 

  /* Set RateData parameters. */
  // NOTE: calc_rates expects these to be the same size as CoolData
  //   RateData.NumberOfTemperatureBins = 600;
  //   RateData.TemperatureStart        = 1;
  //   RateData.TemperatureEnd          = 1.0e8;
  RateData.NumberOfTemperatureBins = CoolData.NumberOfTemperatureBins;
  RateData.TemperatureStart        = CoolData.TemperatureStart;
  RateData.TemperatureEnd          = CoolData.TemperatureEnd;

   
  /* Allocate space in RateData for rates. */
 
  RateData.k1 = new float[RateData.NumberOfTemperatureBins];
  RateData.k2 = new float[RateData.NumberOfTemperatureBins];
  RateData.k3 = new float[RateData.NumberOfTemperatureBins];
  RateData.k4 = new float[RateData.NumberOfTemperatureBins];
  RateData.k5 = new float[RateData.NumberOfTemperatureBins];
  RateData.k6 = new float[RateData.NumberOfTemperatureBins];
  RateData.k7 = new float[RateData.NumberOfTemperatureBins];
  RateData.k8 = new float[RateData.NumberOfTemperatureBins];
  RateData.k9 = new float[RateData.NumberOfTemperatureBins];
  RateData.k10 = new float[RateData.NumberOfTemperatureBins];
  RateData.k11 = new float[RateData.NumberOfTemperatureBins];
  RateData.k12 = new float[RateData.NumberOfTemperatureBins];
  RateData.k13 = new float[RateData.NumberOfTemperatureBins];
  RateData.k13dd = new float[RateData.NumberOfTemperatureBins*7];
  RateData.k14 = new float[RateData.NumberOfTemperatureBins];
  RateData.k15 = new float[RateData.NumberOfTemperatureBins];
  RateData.k16 = new float[RateData.NumberOfTemperatureBins];
  RateData.k17 = new float[RateData.NumberOfTemperatureBins];
  RateData.k18 = new float[RateData.NumberOfTemperatureBins];
  RateData.k19 = new float[RateData.NumberOfTemperatureBins];
  RateData.k20 = new float[RateData.NumberOfTemperatureBins];
  RateData.k21 = new float[RateData.NumberOfTemperatureBins];
  RateData.k22 = new float[RateData.NumberOfTemperatureBins];
  RateData.k23 = new float[RateData.NumberOfTemperatureBins];
  RateData.k50 = new float[RateData.NumberOfTemperatureBins];
  RateData.k51 = new float[RateData.NumberOfTemperatureBins];
  RateData.k52 = new float[RateData.NumberOfTemperatureBins];
  RateData.k53 = new float[RateData.NumberOfTemperatureBins];
  RateData.k54 = new float[RateData.NumberOfTemperatureBins];
  RateData.k55 = new float[RateData.NumberOfTemperatureBins];
  RateData.k56 = new float[RateData.NumberOfTemperatureBins];
  RateData.h2dust = new float[RateData.NumberOfTemperatureBins * 
			      RateData.NumberOfDustTemperatureBins];
  RateData.n_cr_n = new float[RateData.NumberOfTemperatureBins];
  RateData.n_cr_d1 = new float[RateData.NumberOfTemperatureBins];
  RateData.n_cr_d2 = new float[RateData.NumberOfTemperatureBins]; 

  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }
  float afloat = float(a);
  int ioutput = (( MyProcessorNumber == ROOT_PROCESSOR ) ? 1 : 0);
 
  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates)(
     &CoolData.NumberOfTemperatureBins, &afloat, &CoolData.TemperatureStart,
        &CoolData.TemperatureEnd, &CoolData.alpha0, &CoolData.f3,
        &RadiationFieldType, &RateData.CaseBRecombination, &ThreeBodyRate,
     &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
        CoolData.ciHeI,
     CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII,
        CoolData.reHeII1,
     CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp, 
     &PhotoelectricHeatingRate, &CoolData.gammah,
     &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
     CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, CoolData.roth,
        CoolData.rotl,
     CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit,
        CoolData.HDlte, CoolData.HDlow, CoolData.HDcool, CoolData.cieco,
     CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
        CoolData.GAel, CoolData.gas_grain, 
     RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5,
        RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
     RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14,
        RateData.k15, RateData.k16, RateData.k17, RateData.k18,
     RateData.k19, RateData.k20, RateData.k21, RateData.k22, RateData.k23,
     &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27, &RateData.k28,
        &RateData.k29, &RateData.k30, &RateData.k31,
     RateData.k50, RateData.k51, RateData.k52, RateData.k53, RateData.k54,
        RateData.k55, RateData.k56, 
     &RateData.NumberOfDustTemperatureBins, &RateData.DustTemperatureStart, 
     &RateData.DustTemperatureEnd, RateData.h2dust, 
     RateData.n_cr_n, RateData.n_cr_d1, RateData.n_cr_d2, &ioutput);

  /* If using tabulated J21 values for Lyman-Werner, initialize. */
  if (TabulatedLWBackground) {
    if (InitializeLymanWernerTable() == FAIL) {
      ENZO_FAIL("Error in InitializeLymanWernerTable.\n");
    }
  }

  /* If using tabulated HM12 photo-ionization/heating rates, initialize. */
  if (RadiationFieldType == 15) {
    if (InitializeHM12Photorates() == FAIL) {
      ENZO_FAIL("Error in InitializeHM12Photorates.\n"); 
    }
  }

  /* Initialize Cloudy cooling, even if not being used. */
  /* If not used, this will just initialize some data structues. */
  if (InitializeCloudyCooling(Time) == FAIL) {
    ENZO_FAIL("Error in InitializeCloudyCooling.");
  }

  /* If table exists, read metal cooling rates */
  if (MetalCooling == JHW_METAL_COOLING)
    if (ReadMetalCoolingRates(TemperatureUnits, LengthUnits, aUnits, 
			      DensityUnits, TimeUnits, afloat) == FAIL) {
      ENZO_FAIL("Error in ReadMetalCoolingRates.\n");

    }

  return SUCCESS;
}
