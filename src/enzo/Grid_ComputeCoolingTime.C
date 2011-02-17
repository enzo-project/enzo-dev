/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Elizabeth Harper-Clark, August 2009
/              added in CoolingModel parameter
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the cooling time
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
 
/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */
 
extern int RadiationFieldRecomputeMetalRates;
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
int GadgetCoolingTime(float *d, float *e, float *ge, 
		      float *u, float *v, float *w,
		      float *cooltime,
		      int *in, int *jn, int *kn, int *iexpand, 
		      hydro_method *imethod, int *idual, int *idim,
		      int *is, int *js, int *ks, int *ie, int *je, 
		      int *ke, float *dt, float *aye,
		      float *fh, float *utem, float *uxyz, 
		      float *uaye, float *urho, float *utim,
		      float *gamma);

extern "C" void FORTRAN_NAME(cool_multi_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	   float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
           hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co,
	int *ipiht, int *igammah,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa,
	   float *ciHeIa,
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
	float *comp_xraya, float *comp_temp,
           float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
           float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha,
	   float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradfield, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, float *photogamma,
    int *ih2optical, int *iciecool, float *ciecoa,
 	int *icmbTfloor, int *iClHeat, int *iClMMW,
 	float *clMetNorm, float *clEleFra, int *clGridRank, int *clGridDim,
 	float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
 	int *clDataSize, float *clCooling, float *clHeating, float *clMMW);
 
extern "C" void FORTRAN_NAME(cool_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod, int *idual, int *idim, int *igammah,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	float *dt, float *aye, float *temstart, float *temend,
	   float *fh, float *utem,
	float *eta1, float *eta2, float *gamma, float *coola, float *gammaha);
 
int grid::ComputeCoolingTime(float *cooling_time)
{
 
  /* Return if this doesn't concern us. */

  if (RadiativeCooling == 0) return SUCCESS;
 
  if (RadiativeCooling == 0) // would not know what a cooling time is
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
 
  /* Compute the size of the fields. */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Find Multi-species fields. */

  DeNum = HINum = HIINum = HeINum = HeIINum = HeIIINum = HMNum = H2INum = 
    H2IINum = DINum = DIINum = HDINum = 0;
 
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
 
  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];
 
  /* Compute the cooling time. */
 
  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  }
  float afloat = float(a);
 
  /* Metal cooling codes. */
 
  int MetalNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) == -1)
    MetalNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling) {
    if (MetalNum == -1) {
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
      MetalCooling = FALSE;
      MetalNum = 0;
    }
  }
 
  /* Calculate the rates due to the radiation field. */
 
  if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
    ENZO_FAIL("Error in RadiationFieldCalculateRates.\n");
  }
 
  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */
 
  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection *
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection *
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection *
                           double(LengthUnits) * CellWidth[0][0];
 
  /* Call the appropriate FORTRAN routine to do the work. */

  if (MultiSpecies) {
    FORTRAN_NAME(cool_multi_time)(
       density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
       BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum],
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum],
       cooling_time,
       GridDimension, GridDimension+1, GridDimension+2,
          &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
          &HydroMethod,
       &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent, &MetalCooling, 
       &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2,
          GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &CoolData.ih2co, &CoolData.ipiht, &PhotoelectricHeating,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd,
       &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
          CoolData.ciHeI,
       CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII,
          CoolData.reHeII1,
       CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp, &CoolData.gammah,
       &CoolData.comp_xray, &CoolData.temp_xray,
          &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
       BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
          BaryonField[MetalNum],
       CoolData.hyd01k, CoolData.h2k01, CoolData.vibh,
          CoolData.roth, CoolData.rotl,
       CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit,
          CoolData.HDlte, CoolData.HDlow,
       CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp,
       CoolData.GAel,
          CoolData.metals, &CoolData.NumberOfElectronFracBins, 
          &CoolData.ElectronFracStart, &CoolData.ElectronFracEnd,
       RadiationData.Spectrum[0], &RadiationFieldType,
          &RadiationData.NumberOfFrequencyBins,
          &RadiationFieldRecomputeMetalRates,
       &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
       &RadiativeTransfer, BaryonField[gammaNum], 
       &H2OpticalDepthApproximation, &CIECooling, CoolData.cieco,
       &CloudyCoolingData.CMBTemperatureFloor,
       &CloudyCoolingData.IncludeCloudyHeating, &CloudyCoolingData.IncludeCloudyMMW,
       &CloudyCoolingData.CloudyMetallicityNormalization,
       &CloudyCoolingData.CloudyElectronFractionFactor,
       &CloudyCoolingData.CloudyCoolingGridRank,
       CloudyCoolingData.CloudyCoolingGridDimension,
       CloudyCoolingData.CloudyCoolingGridParameters[0],
       CloudyCoolingData.CloudyCoolingGridParameters[1],
       CloudyCoolingData.CloudyCoolingGridParameters[2],
       CloudyCoolingData.CloudyCoolingGridParameters[3],
       CloudyCoolingData.CloudyCoolingGridParameters[4],
       &CloudyCoolingData.CloudyDataSize,
       CloudyCoolingData.CloudyCooling, CloudyCoolingData.CloudyHeating,
       CloudyCoolingData.CloudyMeanMolecularWeight);
  } else if (GadgetEquilibriumCooling==1) {  
    int result = GadgetCoolingTime
      (
       density,totalenergy,gasenergy,velocity1,
       velocity2,velocity3,
       cooling_time,
       GridDimension,GridDimension+1,
       GridDimension+2, &ComovingCoordinates, &HydroMethod,
       &DualEnergyFormalism, &GridRank,
       GridStartIndex,GridStartIndex+1,GridStartIndex+2,
       GridEndIndex,GridEndIndex+1,GridEndIndex+2,&dtFixed,
       &afloat,&CoolData.HydrogenFractionByMass,
       &TemperatureUnits,&LengthUnits,
       &aUnits,&DensityUnits,&TimeUnits,&Gamma);
    if (result == FAIL )  {

      ENZO_FAIL("Error in GadgetCoolingTime.  Exiting.");
    }
  } else { // if not multispecies or Gadget cooling, must be generic cooling.
    FORTRAN_NAME(cool_time)(
       BaryonField[DensNum], BaryonField[TENum], BaryonField[GENum],
          BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
          cooling_time,
       GridDimension, GridDimension+1, GridDimension+2,
          &CoolData.NumberOfTemperatureBins, &ComovingCoordinates,
          &HydroMethod,
       &DualEnergyFormalism, &GridRank, &PhotoelectricHeating,
       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
          GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd, &CoolData.HydrogenFractionByMass,
          &TemperatureUnits,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
          CoolData.EquilibriumRate, &CoolData.gammah);
  }
 
  return SUCCESS;
}
