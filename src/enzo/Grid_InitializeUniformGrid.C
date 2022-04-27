/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A UNIFORM POOL OF GAS)
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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
 
int grid::InitializeUniformGrid(float UniformDensity,
				float UniformTotalEnergy,
				float UniformInternalEnergy,
				float UniformVelocity[], 
				float UniformBField[],
				float UniformCR)
{
  /* declarations */
 
  int dim, i, j, k, index, size, field, GCM;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, MetalIaNum, B1Num, B2Num, B3Num, PhiNum, CRNum;

  int CINum, CIINum, OINum, OIINum, SiINum, SiIINum, SiIIINum, CHINum, CH2INum, 
    CH3IINum, C2INum, COINum, HCOIINum, OHINum, H2OINum, O2INum;


  int ExtraField[2];

  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || HydroMethod > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if ( UseMHD ) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
    if( HydroMethod == MHD_RK ){
        FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    }
    if (UsePoissonDivergenceCleaning) {
      FieldType[NumberOfBaryonFields++] = Phi_pField;
    }
  }

  if ( CRModel ) {
    CRNum = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = CRDensity;
  }

  if (WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;


  int colorfields = NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (TestProblemData.MultiSpecies) {
    FieldType[DeNum     = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum     = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum    = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum    = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum   = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum  = NumberOfBaryonFields++] = HeIIIDensity;
    if (TestProblemData.MultiSpecies > 1) {
      FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (TestProblemData.MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  //  Metal fields, including the standard 'metallicity' as well 
  // as two extra fields
  if (TestProblemData.UseMetallicityField) {
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

    if (StarMakerTypeIaSNe)
      FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

    if(TestProblemData.MultiMetals){
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
  }
 
  // Simon glover's chemistry models (there are several)
  //
  // model #1:  primordial (H, D, He)
  // model #2:  low Z (H, D, He, C, O, Si)
  // model #3:  molecular (H, D, He, C, O, Si, plus molecules)
  // model #4:  simple (just hydrogen)
  // model #5: slightly less simple ( hydrogen plus CO)
  // model #7:  GMC (H, D,He, C, O, plus molecules)
  if(TestProblemData.GloverChemistryModel){

    GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

    // a few species are in all of the models
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;

    // more primordial species
    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
      FieldType[DINum    = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum   = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum   = NumberOfBaryonFields++] = HDIDensity;
      FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
      FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
      FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    }

    // CO molecule
    if( (GCM==3) || (GCM==5) || (GCM==7) ){
      FieldType[COINum   = NumberOfBaryonFields++] = COIDensity;
    }

    // atomic carbon, oxygen
    if( (GCM==2) || (GCM==3) || (GCM==7) ){
      FieldType[CINum     = NumberOfBaryonFields++] = CIDensity;
      FieldType[CIINum    = NumberOfBaryonFields++] = CIIDensity;
      FieldType[OINum     = NumberOfBaryonFields++] = OIDensity;
      FieldType[OIINum    = NumberOfBaryonFields++] = OIIDensity;
    }

    // atomic silicon
    if( (GCM==2) || (GCM==3) ){
      FieldType[SiINum    = NumberOfBaryonFields++] = SiIDensity;
      FieldType[SiIINum   = NumberOfBaryonFields++] = SiIIDensity;
      FieldType[SiIIINum  = NumberOfBaryonFields++] = SiIIIDensity;
    }

    // a ton of molecules
    if( (GCM==3) || (GCM==7) ){
      FieldType[CHINum   = NumberOfBaryonFields++] = CHIDensity;
      FieldType[CH2INum  = NumberOfBaryonFields++] = CH2IDensity;
      FieldType[CH3IINum = NumberOfBaryonFields++] = CH3IIDensity;
      FieldType[C2INum   = NumberOfBaryonFields++] = C2IDensity;
      FieldType[HCOIINum = NumberOfBaryonFields++] = HCOIIDensity;
      FieldType[OHINum   = NumberOfBaryonFields++] = OHIDensity;
      FieldType[H2OINum  = NumberOfBaryonFields++] = H2OIDensity;
      FieldType[O2INum   = NumberOfBaryonFields++] = O2IDensity;
    }

  } //   if(TestProblemData.GloverChemistryModel)

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* compute size of fields */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* allocate fields */
 
  this->AllocateGrids();

  /* set density, total energy */
 
  for (i = 0; i < size; i++) {
    BaryonField[0][i] = UniformDensity;
    BaryonField[1][i] = UniformTotalEnergy;
    if ( CRModel ) BaryonField[CRNum][i] = UniformCR;
  }
 
  /* set velocities */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = UniformVelocity[dim];
 
  /* Set internal energy if necessary. */
 
  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      BaryonField[2][i] = UniformInternalEnergy;

  if (UseMHD) {
    for (dim = 0; dim < 3; dim++) 
      for (i = 0; i < size; i++)
        BaryonField[B1Num+dim][i] = UniformBField[dim];
    if (HydroMethod == MHD_RK){
      for (i=0; i < size; i++)
        BaryonField[PhiNum][i] = 0.;
    }
  }

   /* set density of color fields to user-specified values (if user doesn't specify, 
     the defaults are set in SetDefaultGlobalValues.  Do some minimal amount of error
     checking to try to ensure charge conservation when appropriate */
  for (i = 0; i < size; i++){

    // Set multispecies fields!
    // this attempts to set them such that species conservation is maintained,
    // using the method in CosmologySimulationInitializeGrid.C
    if(TestProblemData.MultiSpecies){

      BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
	TestProblemData.HydrogenFractionByMass * UniformDensity;
 
      BaryonField[HeIINum][i] =  TestProblemData.HeII_Fraction *
	UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIIINum][i] = TestProblemData.HeIII_Fraction *
	UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeINum][i] =
	(1.0 - TestProblemData.HydrogenFractionByMass)*UniformDensity -
	BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];

      if(TestProblemData.MultiSpecies > 1){
	BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
	  BaryonField[HIINum][i];

	BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
	  BaryonField[0][i] * TestProblemData.HydrogenFractionByMass;

	BaryonField[H2IINum][i] = TestProblemData.H2II_Fraction * 2.0 *
	  BaryonField[HIINum][i];
      }

      // HI density is calculated by subtracting off the various ionized fractions
      // from the total
      BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass*BaryonField[0][i]
	- BaryonField[HIINum][i];
      if (MultiSpecies > 1)
	BaryonField[HINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i]
				  + BaryonField[H2INum][i]);

      // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
      // density for convenience) is calculated by summing up all of the ionized species.
      // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
      // calculating mass density, not number density (because the BaryonField values are 4x as
      // heavy for helium for a single electron)
      BaryonField[DeNum][i] = BaryonField[HIINum][i] +
	0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
	  BaryonField[HMNum][i];

      // Set deuterium species (assumed to be a negligible fraction of the total, so not
      // counted in the conservation)
      if(TestProblemData.MultiSpecies > 2){
	BaryonField[DINum ][i]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
	BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
	BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    } // if(TestProblemData.MultiSpecies)

    // metallicity fields (including 'extra' metal fields)
    if(TestProblemData.UseMetallicityField){
      BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* UniformDensity;

      if (StarMakerTypeIaSNe)
	BaryonField[MetalIaNum][i] = TestProblemData.MetallicitySNIaField_Fraction*
	  UniformDensity;

      if(TestProblemData.MultiMetals){
      BaryonField[ExtraField[0]][i] = TestProblemData.MultiMetalsField1_Fraction* UniformDensity;
      BaryonField[ExtraField[1]][i] = TestProblemData.MultiMetalsField2_Fraction* UniformDensity;

      }
    } // if(TestProblemData.UseMetallicityField)

        // simon glover chemistry stuff
    if(TestProblemData.GloverChemistryModel){
      float tempHM, tempH2II;

      GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

      BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
	TestProblemData.HydrogenFractionByMass * UniformDensity;
      BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
	  BaryonField[0][i] * TestProblemData.HydrogenFractionByMass;

      tempHM = TestProblemData.HM_Fraction * BaryonField[HIINum][i];

      tempH2II = TestProblemData.H2II_Fraction * 2.0 * BaryonField[HIINum][i];

      BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass*BaryonField[0][i]
	- BaryonField[HIINum][i];
      BaryonField[HINum][i] -= (tempHM + tempH2II + BaryonField[H2INum][i]);

      if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	BaryonField[DINum   ][i] = TestProblemData.DI_Fraction * UniformDensity;
	BaryonField[DIINum  ][i] = TestProblemData.DII_Fraction * UniformDensity;
	BaryonField[HDINum  ][i] = TestProblemData.HDI_Fraction * UniformDensity;
	BaryonField[HeIINum][i] =  TestProblemData.HeII_Fraction *
	  UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);
	BaryonField[HeIIINum][i] = TestProblemData.HeIII_Fraction *
	  UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);
	BaryonField[HeINum][i] =
	  (1.0 - TestProblemData.HydrogenFractionByMass)*UniformDensity -
	  BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];
      }

      if( (GCM==3) || (GCM==5) || (GCM==7) )
	BaryonField[COINum  ][i] = TestProblemData.COI_Fraction * UniformDensity;

      if( (GCM==2) || (GCM==3) || (GCM==7) ){
	BaryonField[CINum ][i] = TestProblemData.CI_Fraction * UniformDensity;
	BaryonField[CIINum][i] = TestProblemData.CII_Fraction * UniformDensity;
	BaryonField[OINum ][i] = TestProblemData.OI_Fraction * UniformDensity;
	BaryonField[OIINum][i] = TestProblemData.OII_Fraction * UniformDensity;
      }

      if( (GCM==2) || (GCM==3) ){
	BaryonField[SiINum  ][i] = TestProblemData.SiI_Fraction * UniformDensity;
	BaryonField[SiIINum ][i] = TestProblemData.SiII_Fraction * UniformDensity;
	BaryonField[SiIIINum][i] = TestProblemData.SiIII_Fraction * UniformDensity;
      }

      if( (GCM==3) || (GCM==7) ){
	BaryonField[CHINum  ][i] = TestProblemData.CHI_Fraction * UniformDensity;
	BaryonField[CH2INum ][i] = TestProblemData.CH2I_Fraction * UniformDensity;
	BaryonField[CH3IINum][i] = TestProblemData.CH3II_Fraction * UniformDensity;
	BaryonField[C2INum  ][i] = TestProblemData.C2I_Fraction * UniformDensity;
	BaryonField[HCOIINum][i] = TestProblemData.HCOII_Fraction * UniformDensity;
	BaryonField[OHINum  ][i] = TestProblemData.OHI_Fraction * UniformDensity;
	BaryonField[H2OINum ][i] = TestProblemData.H2OI_Fraction * UniformDensity;
	BaryonField[O2INum  ][i] = TestProblemData.O2I_Fraction * UniformDensity;
      }

    } // if(TestProblemData.GloverChemistryModel)
    

  } // for (i = 0; i < size; i++)


  if(UseMHDCT == TRUE){
    for(field=0;field<3;field++)
      for(k=0; k<MagneticDims[field][2]; k++)
	for(j=0; j<MagneticDims[field][1]; j++)
	  for(i=0; i<MagneticDims[field][0];i++){
	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    MagneticField[field][index] = UniformBField[field];
	  }
  }  // if(UseMHDCT == TRUE)

  return SUCCESS;
}
