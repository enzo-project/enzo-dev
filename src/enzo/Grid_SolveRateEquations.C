/***********************************************************************
/
/  GRID CLASS (SOLVE THE MULTI-SPECIES RATE EQUATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
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
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
extern "C" void FORTRAN_NAME(solve_rate)(
	float *de, float *HI, float *HII, float *HeI, float *HeII,
	   float *HeIII, float *tgas, float *dens,
	int *in, int *jn, int *kn, int *nratec, int *iuse_hm,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	float *dt, float *aye, float *temstart, float *temend,
	float *uaye, float *utim, float *urho, float *uxyz,
           float *fh, float *dtoh,
	float *k1a, float *k2a, float *k3a, float *k4a, float *k5a,
	   float *k6a, float *k7a, float *k8a, float *k9a, float *k10a,
	float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a,
           float *k15a,
        float *k16a, float *k17a, float *k18a, float *k19a, float *k22a,
	float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
	   float *k30, float *k31,
	float *k50a, float *k51a, float *k52a, float *k53a, float *k54a,
	   float *k55a, float *k56a,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
	int *iradshield, float *avgsigh, float *avgsighe, float *avgsighe2,
	int *iradfield, float *piHI, float *piHeI,
	int *iradtrans, int *irt_honly, float *kphHI, float *kphHeI, 
	float *kphHeII, float *kdissH2I);
 
 
int grid::SolveRateEquations()
{
 
  /* Return if this doesn't concern us. */
  /* We should be calling SolveRateAndCoolingEquations if both are true */
  
  if (MultiSpecies && RadiativeCooling && 
      (MetalCooling != CLOUDY_METAL_COOLING)) return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck("SolveRateEquations");
 
  /* Declarations */
 
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  FLOAT a = 1.0, dadt;
 
  /* Find Multi-species fields. */
 
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }
 
  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  /* Find the density field. */
 
  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);

  int MetalNum = 0, MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) == -1)
    MetalNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1);
 
  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }
  float afloat = float(a);
 
  /* Calculate the rates due to the radiation field. */
 
  if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
    ENZO_FAIL("Error in RadiationFieldCalculateRates.\n");
  }
 
  /* Set up information for rates which depend on the radiation field. */
 
  /* Precompute factors for self shielding (this is the cross section * dx).
     The factor ToCGS convert code number densities to real number densites
     (the factor of a^3 comes from cancelling a factor of 1/a^3 in
      multi_cool.src) */
 
  float ToCGS = double(DensityUnits)*a*a*a/double(1.67e-24);
  float HIShieldFactor = RadiationData.HIAveragePhotoionizationCrossSection *
                   double(LengthUnits) * CellWidth[0][0] * ToCGS;
  float HeIShieldFactor = RadiationData.HeIAveragePhotoionizationCrossSection*
                    double(LengthUnits) * CellWidth[0][0] * ToCGS;
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoionizationCrossSection *
                     double(LengthUnits) * CellWidth[0][0] * ToCGS;
 
  /* Allocate space for the temperature and compute it. */
 
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float *temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");

  }
 
  /* Call the fortran routine to solve cooling equations. */
  int addRT = (RadiativeTransfer) || (RadiativeTransferFLD);
 
  FORTRAN_NAME(solve_rate)(
       BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum],
          BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum],
          temperature, BaryonField[DensNum],
       GridDimension, GridDimension+1, GridDimension+2,
          &CoolData.NumberOfTemperatureBins, &MultiSpecies,
       GridStartIndex, GridStartIndex+1, GridStartIndex+2,
          GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &dtFixed, &afloat, &CoolData.TemperatureStart,
          &CoolData.TemperatureEnd,
       &aUnits, &TimeUnits, &DensityUnits, &LengthUnits,
          &CoolData.HydrogenFractionByMass,
          &CoolData.DeuteriumToHydrogenRatio,
       RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5,
          RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
       RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14,
          RateData.k15, RateData.k16,
       RateData.k17, RateData.k18, RateData.k19, RateData.k22,
       &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27,
          &RateData.k28, &RateData.k29, &RateData.k30, &RateData.k31,
       RateData.k50, RateData.k51, RateData.k52, RateData.k53,
          RateData.k54, RateData.k55, RateData.k56,
       BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
          BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
       &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
       &RadiationFieldType, &CoolData.piHI, &CoolData.piHeI,
       &addRT, &RadiativeTransferHydrogenOnly,
       BaryonField[kphHINum], BaryonField[kphHeINum], 
       BaryonField[kphHeIINum], BaryonField[kdissH2INum]);
 
  /* deallocate temporary space for solver */
 
  delete temperature;
 
  return SUCCESS;
 
}
