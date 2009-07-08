/***********************************************************************
/
/  INITIALIZE THE MULTI-SPECIES RATES
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:
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
int ReadMetalCoolingRates(float TemperatureUnits, float LengthUnits, 
			  float aUnits, float DensityUnits, float TimeUnits, 
			  float aye);
extern "C" void FORTRAN_NAME(calc_rates)(
     int *nratec, float *aye, float *temstart, float *temend, float *alpha0,
        float *f3, int *iradtype,
     float *utem, float *uxyz, float *uaye, float *urho, float *utim,
     float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
     float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
     float *reHeII2a, float *reHeIIIa, float *brema, float *compa,
     float *piHI, float *piHeI, float *piHeII,
     float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
     float *gpldl, float *gphdl, float *hdlte, float *hdlow,
     float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
     float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, float *k6a,
        float *k7a, float *k8a, float *k9a, float *k10a,
     float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a,
        float *k15a, float *k16a, float *k17a, float *k18a,
     float *k19a, float *k20a, float *k21a, float *k22,
     float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
        float *k30, float *k31,
     float *k50, float *k51, float *k52, float *k53, float *k54, float *k55,
        float *k56, int *ioutput);


// character strings
EXTERN char outfilename[];

 
int InitializeRateData(FLOAT Time)
{
 
  /* Declarations. */
 
  FLOAT a = 1, dadt;
 
  /* Set radiation parameters (ignored if RadiationFieldType = 10 or 11).
     Note: defaults are over-written by CoolData.ParameterFile entries. */
 
  // NOTE: these two are already in the main parameter file
  //  CoolData.alpha0                   = 1.5;     // spectral slope
  //  CoolData.f3                       = 1.0e-21; // spectrum normalization
  CoolData.f0to3                    = 0.1;
  CoolData.RadiationRedshiftOn      = 7.0;
  CoolData.RadiationRedshiftOff     = 0.0;
  CoolData.RadiationRedshiftFullOn  = 6.0;
  CoolData.RadiationRedshiftDropOff = 0.0;
  CoolData.HydrogenFractionByMass   = 0.76;
  /* The DToHRatio is by mass in the code, so multiply by 2. */
  CoolData.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998
  CoolData.NumberOfTemperatureBins = 600;
  CoolData.ih2co                   = 1;
  CoolData.ipiht                   = 1;
  CoolData.TemperatureStart        = 1.0;
  CoolData.TemperatureEnd          = 1.0e8;
  CoolData.comp_xray               = 0;
  CoolData.temp_xray               = 0;


  // over-write defaults if CoolData Parameter file is provided
  if (CoolData.ParameterFilename != NULL) {

    // open parameter file
    FILE *fptr;
    if ((fptr = fopen(CoolData.ParameterFilename, "r")) == NULL) 
      fprintf(stderr,"Error opening CoolData file, keeping defaults.\n");

    else {

      char line[MAX_LINE_LENGTH];
      int ret;

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "RadiationSpectrumNormalization = %"FSYM,
		      &CoolData.f3);
	ret += sscanf(line, "RadiationSpectrumSlope = %"FSYM,
		      &CoolData.alpha0);
	ret += sscanf(line, "CoolDataf0to3 = %"FSYM,
		      &CoolData.f0to3);
	ret += sscanf(line, "RadiationRedshiftOn = %"FSYM,
		      &CoolData.RadiationRedshiftOn);
	ret += sscanf(line, "RadiationRedshiftOff = %"FSYM,
		      &CoolData.RadiationRedshiftOff);
	ret += sscanf(line, "RadiationRedshiftFullOn = %"FSYM,
		      &CoolData.RadiationRedshiftFullOn);
	ret += sscanf(line, "RadiationRedshiftDropOff = %"FSYM,
		      &CoolData.RadiationRedshiftDropOff);
	ret += sscanf(line, "HydrogenFractionByMass = %"FSYM,
		      &CoolData.HydrogenFractionByMass);
	ret += sscanf(line, "DeuteriumToHydrogenRatio = %"FSYM,
		      &CoolData.DeuteriumToHydrogenRatio);
	ret += sscanf(line, "NumberOfTemperatureBins = %"ISYM,
		      &CoolData.NumberOfTemperatureBins);
	ret += sscanf(line, "CoolDataIh2co = %"ISYM, &CoolData.ih2co);
	ret += sscanf(line, "CoolDataIpiht = %"ISYM, &CoolData.ipiht);
	ret += sscanf(line, "TemperatureStart = %"FSYM,
		      &CoolData.TemperatureStart);
	ret += sscanf(line, "TemperatureEnd = %"FSYM,
		      &CoolData.TemperatureEnd);
	ret += sscanf(line, "CoolDataCompXray = %"FSYM,
		      &CoolData.comp_xray);
	ret += sscanf(line, "CoolDataTempXray = %"FSYM,
		      &CoolData.temp_xray);
      }

      // clean up
      rewind(fptr);

      // close parameter file
      fclose(fptr);
    }
  }  // end reading cooldata parameter file 
 
  if (debug) printf("InitializeRateData: NumberOfTemperatureBins = %"ISYM"\n",
		    CoolData.NumberOfTemperatureBins);

  
 
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
  CoolData.GAHI    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAH2    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAHe    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAHp    = new float[CoolData.NumberOfTemperatureBins];
  CoolData.GAel    = new float[CoolData.NumberOfTemperatureBins];
 
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
  RateData.k50 = new float[RateData.NumberOfTemperatureBins];
  RateData.k51 = new float[RateData.NumberOfTemperatureBins];
  RateData.k52 = new float[RateData.NumberOfTemperatureBins];
  RateData.k53 = new float[RateData.NumberOfTemperatureBins];
  RateData.k54 = new float[RateData.NumberOfTemperatureBins];
  RateData.k55 = new float[RateData.NumberOfTemperatureBins];
  RateData.k56 = new float[RateData.NumberOfTemperatureBins];
 
  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      ENZO_FAIL("");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }
  float afloat = float(a);
  int ioutput = (( MyProcessorNumber == ROOT_PROCESSOR ) ? 1 : 0);
 
  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates)(
     &CoolData.NumberOfTemperatureBins, &afloat, &CoolData.TemperatureStart,
        &CoolData.TemperatureEnd, &CoolData.alpha0, &CoolData.f3,
        &RadiationFieldType,
     &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
        CoolData.ciHeI,
     CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII,
        CoolData.reHeII1,
     CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp,
     &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
     CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, CoolData.roth,
        CoolData.rotl,
     CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit,
        CoolData.HDlte, CoolData.HDlow,
     CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
        CoolData.GAel,
     RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5,
        RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
     RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14,
        RateData.k15, RateData.k16, RateData.k17, RateData.k18,
     RateData.k19, RateData.k20, RateData.k21, RateData.k22,
     &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27, &RateData.k28,
        &RateData.k29, &RateData.k30, &RateData.k31,
     RateData.k50, RateData.k51, RateData.k52, RateData.k53, RateData.k54,
        RateData.k55, RateData.k56, &ioutput);
 

  // output CoolData parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s.\n", 
	      outfilename);
      ENZO_FAIL("");
    }

    fprintf(outfptr,"RadiationSpectrumNormalization = %"FSYM"\n",
	    CoolData.f3);
    fprintf(outfptr,"RadiationSpectrumSlope = %"FSYM"\n",
	    CoolData.alpha0);
    fprintf(outfptr,"CoolDataf0to3 = %"FSYM"\n", CoolData.f0to3);
    fprintf(outfptr,"RadiationRedshiftOn = %"FSYM"\n",
	    CoolData.RadiationRedshiftOn);
    fprintf(outfptr,"RadiationRedshiftOff = %"FSYM"\n",
	    CoolData.RadiationRedshiftOff);
    fprintf(outfptr,"RadiationRedshiftFullOn = %"FSYM"\n",
	    CoolData.RadiationRedshiftFullOn);
    fprintf(outfptr,"RadiationRedshiftDropOff = %"FSYM"\n",
	    CoolData.RadiationRedshiftDropOff);
    fprintf(outfptr,"HydrogenFractionByMass = %"FSYM"\n",
	    CoolData.HydrogenFractionByMass);
    fprintf(outfptr,"DeuteriumToHydrogenRatio = %"FSYM"\n",
	    CoolData.DeuteriumToHydrogenRatio);
    fprintf(outfptr,"NumberOfTemperatureBins = %"ISYM"\n",
	    CoolData.NumberOfTemperatureBins);
    fprintf(outfptr,"CoolDataIh2co = %"ISYM"\n", CoolData.ih2co);
    fprintf(outfptr,"CoolDataIpiht = %"ISYM"\n", CoolData.ipiht);
    fprintf(outfptr,"TemperatureStart = %"FSYM"\n",
	    CoolData.TemperatureStart);
    fprintf(outfptr,"TemperatureEnd = %"FSYM"\n",
	    CoolData.TemperatureEnd);
    fprintf(outfptr,"CoolDataCompXray = %"FSYM"\n", CoolData.comp_xray);
    fprintf(outfptr,"CoolDataTempXray = %"FSYM"\n", CoolData.temp_xray);
    
    // close parameter file
    fclose(outfptr);
  }

  /* If table exists, read metal cooling rates */

  if (MetalCooling == JHW_METAL_COOLING)
    if (ReadMetalCoolingRates(TemperatureUnits, LengthUnits, aUnits, 
			      DensityUnits, TimeUnits, afloat) == FAIL) {
      fprintf(stderr, "Error in ReadMetalCoolingRates.\n");
      ENZO_FAIL("");
    }

  return SUCCESS;
}
