/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  July, 2005 to solve cool and rate equations simultaneously
/  modified2:  2007-2008 to use new solver
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

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

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
double ReturnWallTime();
extern "C" void FORTRAN_NAME(primordial_solver)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand, 
           hydro_method *imethod,
        int *idual, int *ispecies, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co, 
	   int *ipiht,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim, float *uvel,
	float *eta1, float *eta2, float *gamma, float *fh, float *dtoh,
	float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, 
	   float *k6a, float *k7a, float *k8a, float *k9a, float *k10a,
	float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a, 
           float *k15a,
        float *k16a, float *k17a, float *k18a, float *k19a, float *k21a,
            float *k22a, float *k23a,
	float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
	   float *k30, float *k31,
	float *k50a, float *k51a, float *k52a, float *k53a, float *k54a,
	   float *k55a, float *k56a,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, 
	   float *ciHeIa, 
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a, 
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa,
	float *comp_xraya, float *comp_temp, 
           float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa, float *HDcool, float *ciecoa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *inutot, int *iradtype, int *nfreq, 
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
    int *iciecool, int *ih2optical, int *errcode, int *omaskflag, int *threebody, int *subgridmask
#ifdef UNUSED_TABULATED_EQ
    ,float *HIeqtable, float *HIIeqtable, float *H2Ieqtable, 
        int *nrhobins, int *nebins, float *rhostart, float *rhoend,
    float *estart, float *eend
#endif
    );

int grid::SolveHighDensityPrimordialChemistry()
{
  /* Return if this doesn't concern us. */
  if (!( (MultiSpecies == 3)
      &&  RadiativeCooling 
      && (PrimordialChemistrySolver == 1))) return SUCCESS;

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("SolveHighDensityPrimordialChemistry");

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  FLOAT a = 1.0, dadt;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  /* Compute size of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* Get easy to handle pointers for each variable. */

  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  /* Compute total gas energy if using MHD */
  if (HydroMethod == MHD_RK) {
    totalenergy = new float[size];
    float B2;
    for (int n=0; n<size; n++) {
      B2 = pow(BaryonField[B1Num][n],2) + pow(BaryonField[B2Num][n],2) + pow(BaryonField[B3Num][n],2);
      totalenergy[n] = BaryonField[TENum][n] - 0.5*B2/BaryonField[DensNum][n];
    }
  }
  else {
    totalenergy = BaryonField[TENum];
  }


  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, MassUnits = 1, aUnits = 1;

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
            ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  float afloat = float(a);

  /* Calculate the rates due to the radiation field. */

  if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
        ENZO_FAIL("Error in RadiationFieldCalculateRates.");
  }

  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */

  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection * 
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection * 
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection * 
                           double(LengthUnits) * CellWidth[0][0];

  /* Call the fortran routine to solve cooling equations. */

  int RTCoupledSolverIntermediateStep = FALSE;
  int mask = 1;
  int ErrCode = 0, OutputFlag = 0;

  FORTRAN_NAME(primordial_solver)(
    density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
    BaryonField[DeNum], BaryonField[HINum], BaryonField[HIINum], 
       BaryonField[HeINum], BaryonField[HeIINum], BaryonField[HeIIINum], 
    GridDimension, GridDimension+1, GridDimension+2, 
       &CoolData.NumberOfTemperatureBins, &ComovingCoordinates, &HydroMethod, 
    &DualEnergyFormalism, &MultiSpecies, &GridRank,
    GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
       GridEndIndex, GridEndIndex+1, GridEndIndex+2,
       &CoolData.ih2co, &CoolData.ipiht,
    &dtFixed, &afloat, &CoolData.TemperatureStart, &CoolData.TemperatureEnd,
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits, &VelocityUnits,
    &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
       &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio,
    RateData.k1, RateData.k2, RateData.k3, RateData.k4, RateData.k5, 
       RateData.k6, RateData.k7, RateData.k8, RateData.k9, RateData.k10,
    RateData.k11, RateData.k12, RateData.k13, RateData.k13dd, RateData.k14, 
       RateData.k15, RateData.k16,
    RateData.k17, RateData.k18, RateData.k19, RateData.k21, RateData.k22,
       RateData.k23,
    &RateData.k24, &RateData.k25, &RateData.k26, &RateData.k27,
       &RateData.k28, &RateData.k29, &RateData.k30, &RateData.k31,
    RateData.k50, RateData.k51, RateData.k52, RateData.k53,
       RateData.k54, RateData.k55, RateData.k56,
    CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI,
       CoolData.ciHeI, 
    CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, 
    CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, &CoolData.comp,
    &CoolData.comp_xray, &CoolData.temp_xray,
       &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
    BaryonField[HMNum], BaryonField[H2INum], BaryonField[H2IINum],
       BaryonField[DINum], BaryonField[DIINum], BaryonField[HDINum],
    CoolData.hyd01k, CoolData.h2k01, CoolData.vibh, CoolData.roth,CoolData.rotl,
    CoolData.GP99LowDensityLimit, CoolData.GP99HighDensityLimit, 
       CoolData.HDlte, CoolData.HDlow, CoolData.HDcool, CoolData.cieco,
    CoolData.GAHI, CoolData.GAH2, CoolData.GAHe, CoolData.GAHp, CoolData.GAel,
    RadiationData.Spectrum[0], &RadiationFieldType, 
          &RadiationData.NumberOfFrequencyBins, 
    &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
    &CIECooling, &H2OpticalDepthApproximation, &ErrCode, &OutputFlag,
           &ThreeBodyRate, &mask
#ifdef UNUSED_TABULATED_EQ
        ,RateData.HighDensityEquilibriumRate[0], 
        RateData.HighDensityEquilibriumRate[1], 
        RateData.HighDensityEquilibriumRate[2], 
    &RateData.HighDensityNumberOfDensityBins, &RateData.HighDensityNumberOfEnergyBins,
    &RateData.HighDensityDensityStart, &RateData.HighDensityDensityStop,
    &RateData.HighDensityEnergyStart, &RateData.HighDensityEnergyStop
#endif
    );

  if (ErrCode) {
      fprintf(stdout, "GridLeftEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(stdout, "GridRightEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
      fprintf(stdout, "GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	      GridDimension[0], GridDimension[1], GridDimension[2]);
      ENZO_FAIL("Error in FORTRAN rate/cool solver!\n");
  }

  if (HydroMethod == MHD_RK) {
    float B2, v2;
    for (int n = 0; n < size; n++) {
      B2 = pow(BaryonField[B1Num][n],2) + pow(BaryonField[B2Num][n],2) + pow(BaryonField[B3Num][n],2);

      /* Always trust gas energy in cooling routine */
      if (DualEnergyFormalism) {

	v2 = pow(BaryonField[Vel1Num][n],2) + 
	  pow(BaryonField[Vel2Num][n],2) + pow(BaryonField[Vel3Num][n],2);
	BaryonField[TENum][n] = gasenergy[n] + 0.5*v2 + 0.5*B2/BaryonField[DensNum][n];
      }
      else {
	BaryonField[TENum][n] = totalenergy[n] + 0.5*B2/BaryonField[DensNum][n];
      }
      
    }
    
    delete totalenergy;
  }

  return SUCCESS;

}
