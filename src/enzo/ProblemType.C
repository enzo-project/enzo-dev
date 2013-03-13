/***********************************************************************
/
/  PROBLEM TYPE CLASS
/
/  written by: Matthew Turk, Oliver Hahn
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#include "ProblemType.h"

EnzoProblemMap& get_problem_types()
{
    static EnzoProblemMap problem_type_map;
    return problem_type_map;
}

/* This takes a string, grabs the (static) plugin map defined above, and
   returns the plugin creator for that. */
EnzoProblemType *select_problem_type( std::string problem_type_name)
{
    EnzoProblemType_creator *ept_creator = get_problem_types()
            [ problem_type_name ];

    /* Simply throw an error if no such plugin exists... */

    if( !ept_creator )
    {   
        EnzoProblemMap mymap = get_problem_types();

        for (EnzoProblemMap::const_iterator it 
                = mymap.begin(); it != mymap.end(); ++it) {
          std::cout << "Available: " << it->first << std::endl;
        }
        ENZO_FAIL("Unknown output plug-in.");
    }

    EnzoProblemType *ptype = ept_creator->create();

    return ptype;

}

EnzoProblemType::EnzoProblemType()
{
    this->DataLabelCount = 0;
}

int EnzoProblemType::AddDataLabel(const char *FieldName) {
    /* We allocate a new copy of FieldName */
    /* Include NUL-terminator */
    int slen = strlen(FieldName) + 1;
    char *fcopy = new char[slen];
    strncpy(fcopy, FieldName, slen);
    DataLabel[this->DataLabelCount] = fcopy;
    std::cout << "Adding " << DataLabel[this->DataLabelCount]
              << " in place " << this->DataLabelCount << std::endl;
    return DataLabelCount++;
}

grid *EnzoProblemType::CreateNewUniformGrid(
                grid *ParentGrid,
                int Rank, int Dimensions[],
                FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles,
                float UniformDensity,
				float UniformTotalEnergy,
				float UniformInternalEnergy,
				float UniformVelocity[], 
				float UniformBField[])
{
    grid *grid_data = new grid;
    grid_data->InheritProperties(ParentGrid);
    grid_data->PrepareGrid(Rank, Dimensions, LeftEdge, RightEdge,
            NumParticles);
    this->InitializeUniformGrid(grid_data,
                UniformDensity,
				UniformTotalEnergy,
				UniformInternalEnergy,
				UniformVelocity, 
				UniformBField);
    return grid_data;
}

int EnzoProblemType::InitializeUniformGrid(
                grid *tg,
                float UniformDensity,
				float UniformTotalEnergy,
				float UniformInternalEnergy,
				float UniformVelocity[], 
				float UniformBField[])
{
  /* declarations */
 
  int dim, i, size, field, GCM;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, B1Num, B2Num, B3Num, PhiNum;

  int CINum, CIINum, OINum, OIINum, SiINum, SiIINum, SiIIINum, CHINum, CH2INum, 
    CH3IINum, C2INum, COINum, HCOIINum, OHINum, H2OINum, O2INum;


  int ExtraField[2];

  /* create fields */
 
  tg->NumberOfBaryonFields = 0;
  tg->FieldType[tg->NumberOfBaryonFields++] = Density;
  tg->FieldType[tg->NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    tg->FieldType[tg->NumberOfBaryonFields++] = InternalEnergy;
  int vel = tg->NumberOfBaryonFields;
  tg->FieldType[tg->NumberOfBaryonFields++] = Velocity1;
  if (tg->GridRank > 1 || HydroMethod > 2)
    tg->FieldType[tg->NumberOfBaryonFields++] = Velocity2;
  if (tg->GridRank > 2 || HydroMethod > 2)
    tg->FieldType[tg->NumberOfBaryonFields++] = Velocity3;
  if (HydroMethod == MHD_RK) {
    tg->FieldType[B1Num = tg->NumberOfBaryonFields++] = Bfield1;
    tg->FieldType[B2Num = tg->NumberOfBaryonFields++] = Bfield2;
    tg->FieldType[B3Num = tg->NumberOfBaryonFields++] = Bfield3;
    tg->FieldType[PhiNum = tg->NumberOfBaryonFields++] = PhiField;
    if (UseDivergenceCleaning) {
      tg->FieldType[tg->NumberOfBaryonFields++] = Phi_pField;
    }
  }


  int colorfields = tg->NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (TestProblemData.MultiSpecies) {
    tg->FieldType[DeNum     = tg->NumberOfBaryonFields++] = ElectronDensity;
    tg->FieldType[HINum     = tg->NumberOfBaryonFields++] = HIDensity;
    tg->FieldType[HIINum    = tg->NumberOfBaryonFields++] = HIIDensity;
    tg->FieldType[HeINum    = tg->NumberOfBaryonFields++] = HeIDensity;
    tg->FieldType[HeIINum   = tg->NumberOfBaryonFields++] = HeIIDensity;
    tg->FieldType[HeIIINum  = tg->NumberOfBaryonFields++] = HeIIIDensity;
    if (TestProblemData.MultiSpecies > 1) {
      tg->FieldType[HMNum   = tg->NumberOfBaryonFields++] = HMDensity;
      tg->FieldType[H2INum  = tg->NumberOfBaryonFields++] = H2IDensity;
      tg->FieldType[H2IINum = tg->NumberOfBaryonFields++] = H2IIDensity;
    }
    if (TestProblemData.MultiSpecies > 2) {
      tg->FieldType[DINum   = tg->NumberOfBaryonFields++] = DIDensity;
      tg->FieldType[DIINum  = tg->NumberOfBaryonFields++] = DIIDensity;
      tg->FieldType[HDINum  = tg->NumberOfBaryonFields++] = HDIDensity;
    }
  }

  //  Metal fields, including the standard 'metallicity' as well 
  // as two extra fields
  if (TestProblemData.UseMetallicityField) {
    tg->FieldType[MetalNum = tg->NumberOfBaryonFields++] = Metallicity;

    if(TestProblemData.MultiMetals){
      tg->FieldType[ExtraField[0] = tg->NumberOfBaryonFields++] = ExtraType0;
      tg->FieldType[ExtraField[1] = tg->NumberOfBaryonFields++] = ExtraType1;
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
    tg->FieldType[HIINum   = tg->NumberOfBaryonFields++] = HIIDensity;
    tg->FieldType[HINum    = tg->NumberOfBaryonFields++] = HIDensity;
    tg->FieldType[H2INum   = tg->NumberOfBaryonFields++] = H2IDensity;

    // more primordial species
    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
      tg->FieldType[DINum    = tg->NumberOfBaryonFields++] = DIDensity;
      tg->FieldType[DIINum   = tg->NumberOfBaryonFields++] = DIIDensity;
      tg->FieldType[HDINum   = tg->NumberOfBaryonFields++] = HDIDensity;
      tg->FieldType[HeINum   = tg->NumberOfBaryonFields++] = HeIDensity;
      tg->FieldType[HeIINum  = tg->NumberOfBaryonFields++] = HeIIDensity;
      tg->FieldType[HeIIINum = tg->NumberOfBaryonFields++] = HeIIIDensity;
    }

    // CO molecule
    if( (GCM==3) || (GCM==5) || (GCM==7) ){
      tg->FieldType[COINum   = tg->NumberOfBaryonFields++] = COIDensity;
    }

    // atomic carbon, oxygen
    if( (GCM==2) || (GCM==3) || (GCM==7) ){
      tg->FieldType[CINum     = tg->NumberOfBaryonFields++] = CIDensity;
      tg->FieldType[CIINum    = tg->NumberOfBaryonFields++] = CIIDensity;
      tg->FieldType[OINum     = tg->NumberOfBaryonFields++] = OIDensity;
      tg->FieldType[OIINum    = tg->NumberOfBaryonFields++] = OIIDensity;
    }

    // atomic silicon
    if( (GCM==2) || (GCM==3) ){
      tg->FieldType[SiINum    = tg->NumberOfBaryonFields++] = SiIDensity;
      tg->FieldType[SiIINum   = tg->NumberOfBaryonFields++] = SiIIDensity;
      tg->FieldType[SiIIINum  = tg->NumberOfBaryonFields++] = SiIIIDensity;
    }

    // a ton of molecules
    if( (GCM==3) || (GCM==7) ){
      tg->FieldType[CHINum   = tg->NumberOfBaryonFields++] = CHIDensity;
      tg->FieldType[CH2INum  = tg->NumberOfBaryonFields++] = CH2IDensity;
      tg->FieldType[CH3IINum = tg->NumberOfBaryonFields++] = CH3IIDensity;
      tg->FieldType[C2INum   = tg->NumberOfBaryonFields++] = C2IDensity;
      tg->FieldType[HCOIINum = tg->NumberOfBaryonFields++] = HCOIIDensity;
      tg->FieldType[OHINum   = tg->NumberOfBaryonFields++] = OHIDensity;
      tg->FieldType[H2OINum  = tg->NumberOfBaryonFields++] = H2OIDensity;
      tg->FieldType[O2INum   = tg->NumberOfBaryonFields++] = O2IDensity;
    }

  } //   if(TestProblemData.GloverChemistryModel)

  /* Return if this doesn't concern us. */
 
  if (tg->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* compute size of fields */
 
  size = 1;
  for (dim = 0; dim < tg->GridRank; dim++)
    size *= tg->GridDimension[dim];
 
  /* allocate fields */
 
  // for (field = 0; field < tg->NumberOfBaryonFields; field++)
  //   if (tg->BaryonField[field] == NULL)
  //     tg->BaryonField[field] = new float[size];
  tg->AllocateGrids();  
  /* set density, total energy */
 
  for (i = 0; i < size; i++) {
    tg->BaryonField[0][i] = UniformDensity;
    tg->BaryonField[1][i] = UniformTotalEnergy;
  }
 
  /* set velocities */
 
  for (dim = 0; dim < tg->GridRank; dim++)
    for (i = 0; i < size; i++)
      tg->BaryonField[vel+dim][i] = UniformVelocity[dim];
 
  /* Set internal energy if necessary. */
 
  if (DualEnergyFormalism)
    for (i = 0; i < size; i++)
      tg->BaryonField[2][i] = UniformInternalEnergy;

  if (HydroMethod == MHD_RK) {
    for (dim = 0; dim < 3; dim++) 
      for (i = 0; i < size; i++)
	tg->BaryonField[B1Num+dim][i] = UniformBField[dim];
  }

   /* set density of color fields to user-specified values (if user doesn't specify, 
     the defaults are set in SetDefaultGlobalValues.  Do some minimal amount of error
     checking to try to ensure charge conservation when appropriate */
  for (i = 0; i < size; i++){

    // Set multispecies fields!
    // this attempts to set them such that species conservation is maintained,
    // using the method in CosmologySimulationInitializeGrid.C
    if(TestProblemData.MultiSpecies){

      tg->BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
	TestProblemData.HydrogenFractionByMass * UniformDensity;
 
      tg->BaryonField[HeIINum][i] =  TestProblemData.HeII_Fraction *
	UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      tg->BaryonField[HeIIINum][i] = TestProblemData.HeIII_Fraction *
	UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);

      tg->BaryonField[HeINum][i] =
	(1.0 - TestProblemData.HydrogenFractionByMass)*UniformDensity -
	tg->BaryonField[HeIINum][i] - tg->BaryonField[HeIIINum][i];

      if(TestProblemData.MultiSpecies > 1){
	tg->BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
	  tg->BaryonField[HIINum][i];

	tg->BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
	  tg->BaryonField[0][i] * TestProblemData.HydrogenFractionByMass;

	tg->BaryonField[H2IINum][i] = TestProblemData.H2II_Fraction * 2.0 *
	  tg->BaryonField[HIINum][i];
      }

      // HI density is calculated by subtracting off the various ionized fractions
      // from the total
      tg->BaryonField[HINum][i] = 
            TestProblemData.HydrogenFractionByMass*tg->BaryonField[0][i]
            - tg->BaryonField[HIINum][i];
      if (MultiSpecies > 1)
	tg->BaryonField[HINum][i] -= (tg->BaryonField[HMNum][i]
                                + tg->BaryonField[H2IINum][i]
				                + tg->BaryonField[H2INum][i]);

      // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
      // density for convenience) is calculated by summing up all of the ionized species.
      // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
      // calculating mass density, not number density (because the BaryonField values are 4x as
      // heavy for helium for a single electron)
      tg->BaryonField[DeNum][i] = tg->BaryonField[HIINum][i] +
	0.25*tg->BaryonField[HeIINum][i] + 0.5*tg->BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	tg->BaryonField[DeNum][i] += 0.5*tg->BaryonField[H2IINum][i] -
	  tg->BaryonField[HMNum][i];

      // Set deuterium species (assumed to be a negligible fraction of the total, so not
      // counted in the conservation)
      if(TestProblemData.MultiSpecies > 2){
	tg->BaryonField[DINum ][i]  = TestProblemData.DeuteriumToHydrogenRatio * tg->BaryonField[HINum][i];
	tg->BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * tg->BaryonField[HIINum][i];
	tg->BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * tg->BaryonField[H2INum][i];
      }

    } // if(TestProblemData.MultiSpecies)

    // metallicity fields (including 'extra' metal fields)
    if(TestProblemData.UseMetallicityField){
      tg->BaryonField[MetalNum][i] = TestProblemData.MetallicityField_Fraction* UniformDensity;

      if(TestProblemData.MultiMetals){
      tg->BaryonField[ExtraField[0]][i] = TestProblemData.MultiMetalsField1_Fraction* UniformDensity;
      tg->BaryonField[ExtraField[1]][i] = TestProblemData.MultiMetalsField2_Fraction* UniformDensity;

      }
    } // if(TestProblemData.UseMetallicityField)

        // simon glover chemistry stuff
    if(TestProblemData.GloverChemistryModel){
      float tempHM, tempH2II;

      GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

      tg->BaryonField[HIINum][i] = TestProblemData.HII_Fraction * 
	TestProblemData.HydrogenFractionByMass * UniformDensity;
      tg->BaryonField[H2INum][i] = TestProblemData.H2I_Fraction *
	  tg->BaryonField[0][i] * TestProblemData.HydrogenFractionByMass;

      tempHM = TestProblemData.HM_Fraction * tg->BaryonField[HIINum][i];

      tempH2II = TestProblemData.H2II_Fraction * 2.0 * tg->BaryonField[HIINum][i];

      tg->BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass*tg->BaryonField[0][i]
	- tg->BaryonField[HIINum][i];
      tg->BaryonField[HINum][i] -= (tempHM + tempH2II + tg->BaryonField[H2INum][i]);

      if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	tg->BaryonField[DINum   ][i] = TestProblemData.DI_Fraction * UniformDensity;
	tg->BaryonField[DIINum  ][i] = TestProblemData.DII_Fraction * UniformDensity;
	tg->BaryonField[HDINum  ][i] = TestProblemData.HDI_Fraction * UniformDensity;
	tg->BaryonField[HeIINum][i] =  TestProblemData.HeII_Fraction *
	  UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);
	tg->BaryonField[HeIIINum][i] = TestProblemData.HeIII_Fraction *
	  UniformDensity * (1.0-TestProblemData.HydrogenFractionByMass);
	tg->BaryonField[HeINum][i] =
	  (1.0 - TestProblemData.HydrogenFractionByMass)*UniformDensity -
	  tg->BaryonField[HeIINum][i] - tg->BaryonField[HeIIINum][i];
      }

      if( (GCM==3) || (GCM==5) || (GCM==7) )
	tg->BaryonField[COINum  ][i] = TestProblemData.COI_Fraction * UniformDensity;

      if( (GCM==2) || (GCM==3) || (GCM==7) ){
	tg->BaryonField[CINum ][i] = TestProblemData.CI_Fraction * UniformDensity;
	tg->BaryonField[CIINum][i] = TestProblemData.CII_Fraction * UniformDensity;
	tg->BaryonField[OINum ][i] = TestProblemData.OI_Fraction * UniformDensity;
	tg->BaryonField[OIINum][i] = TestProblemData.OII_Fraction * UniformDensity;
      }

      if( (GCM==2) || (GCM==3) ){
	tg->BaryonField[SiINum  ][i] = TestProblemData.SiI_Fraction * UniformDensity;
	tg->BaryonField[SiIINum ][i] = TestProblemData.SiII_Fraction * UniformDensity;
	tg->BaryonField[SiIIINum][i] = TestProblemData.SiIII_Fraction * UniformDensity;
      }

      if( (GCM==3) || (GCM==7) ){
	tg->BaryonField[CHINum  ][i] = TestProblemData.CHI_Fraction * UniformDensity;
	tg->BaryonField[CH2INum ][i] = TestProblemData.CH2I_Fraction * UniformDensity;
	tg->BaryonField[CH3IINum][i] = TestProblemData.CH3II_Fraction * UniformDensity;
	tg->BaryonField[C2INum  ][i] = TestProblemData.C2I_Fraction * UniformDensity;
	tg->BaryonField[HCOIINum][i] = TestProblemData.HCOII_Fraction * UniformDensity;
	tg->BaryonField[OHINum  ][i] = TestProblemData.OHI_Fraction * UniformDensity;
	tg->BaryonField[H2OINum ][i] = TestProblemData.H2OI_Fraction * UniformDensity;
	tg->BaryonField[O2INum  ][i] = TestProblemData.O2I_Fraction * UniformDensity;
      }

    } // if(TestProblemData.GloverChemistryModel)

  } // for (i = 0; i < size; i++)

  return SUCCESS;
}

void EnzoProblemType::FinalizeGrids(HierarchyEntry **RefLevels,
            HierarchyEntry &TopGrid, TopGridData &MetaData)
{

  /* set up subgrids from level 1 to max refinement level -1 */
  
  int lev;
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (RefLevels[lev]->GridData->ProjectSolutionToParentGrid(
          *(RefLevels[lev-1]->GridData))
        == FAIL) {
      ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  /* set up the root grid */

  if (MaximumRefinementLevel > 0) {
    if (RefLevels[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
        == FAIL) {
      ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData) == FAIL) {
      ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
    }

}

#endif
