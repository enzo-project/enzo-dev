/**********************************************************
//
//  INITIALIZE CHEMICAL EVOLUTION MODEL TEST
//
//
//  written by : Andrew Emerick
//  date:        January 2016
//  modified:
//
//  PURPOSE:
//   Simple test of coupling star particles to chemical
//   evolution model. Single star particle of desired mass
//   is placed at desired position and evolved. Test of
//   metal ejection via stellar winds over particle lifetime
//   as well as end of life metal ejection (if applicable)
//
//
//
//   RETURNS: SUCCESS or FAIL
//
**********************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <string.h>
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"
#include "TopGridData.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level);

char* ChemicalSpeciesBaryonFieldLabel(const int &atomic_number, int element_set=1);

void WriteListOfFloats(FILE *fptr, int N, float floats[]);

int IndividualStarProperties_Initialize(TopGridData &MetaData);
int IndividualStarRadiationProperties_Initialize(void);
int InitializeStellarYields(const float &time);
int InitializeDoublePowerDarkMatter(void);

int ChemicalEvolutionTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                                    TopGridData &MetaData)
{
  char *DensName  = "Density";
  char *TEName    = "TotalEnergy";
  char *GEName    = "GasEnergy";
  char *Vel1Name  = "x-velocity";
  char *Vel2Name  = "y-velocity";
  char *Vel3Name  = "z-velocity";
  char *CRName      = "CREnergyDensity";
  char *GravPotentialName = "GravPotential";
  char *MetalName   = "Metal_Density";
  char *MetalIaName = "MetalSNIa_Density";

  char *AGBMetalName    = "AGB_Metal_Density";
  char *PopIIIMetalName = "PopIII_Metal_Density";
  char *PopIIIPISNeMetalName = "PopIII_PISNe_Metal_Density";
  char *WindMetalName = "Intermediate_Wind_Metal_Density";
  char *WindMetalName2 = "Massive_Wind_Metal_Density";
  char *SNIIMetalName = "SNII_Metal_Density";
  char *SNIaMetalName = "SNIa_Metal_Density";
  char *RProcMetalName = "RProcess_Metal_Density";
  char *ExtraMetalName0    = "SNIa_sCh_Metal_Density";
  char *ExtraMetalName1    = "SNIa_SDS_Metal_Density";
  char *ExtraMetalName2    = "SNIa_HeRS_Metal_Density";


  /* Names for Primordial chemistry */
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";

  /* Names for chemical evolution element abundances */
     // these are handled in lookup table in function declared inGrid_IdentifyChemicalTracerSpeciesFields //

  /* declarations */
  char line[MAX_LINE_LENGTH];
  int  dim, ret, level, disk, i; // might not need disk

  /* make sure we are in 3D */

  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do ChemicalEvolutionTest in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default values for parameters */
  float ChemicalEvolutionTestGasDensity     = 1.673E-24,
        ChemicalEvolutionTestGasTemperature = 1.0E4,
        ChemicalEvolutionTestGasMetallicity = tiny_number;
  ChemicalEvolutionTestConcentration = 10.0;
  ChemicalEvolutionTestGasRadius      = 0.25;         // code units

  ChemicalEvolutionTestBackgroundGasDensity = ChemicalEvolutionTestGasDensity;
  ChemicalEvolutionTestBackgroundGasTemperature = ChemicalEvolutionTestGasTemperature;

  int   ChemicalEvolutionTestRefineAtStart   = 1,
        ChemicalEvolutionTestUseMetals       = 1;


  ChemicalEvolutionTestGasDistribution = 0;

  /* MultiSpecies parameters. Ionized, primordial gas */
  TestProblemData.HydrogenFractionByMass   = 0.75;
  TestProblemData.HII_Fraction             = 1.00;
  TestProblemData.HeII_Fraction            = 0.00;
  TestProblemData.HeIII_Fraction           = 1.00;
  TestProblemData.HM_Fraction              = 0.00;
  TestProblemData.H2I_Fraction             = 0.00;
  TestProblemData.H2II_Fraction            = 0.00;
  TestProblemData.DeuteriumToHydrogenRatio = 0.00;

  /* MultiMetals Parameters assigned defaults in SetDefault to tiny_number */

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;

    ret += sscanf(line, "ChemicalEvolutionTestNumberOfStars = %"ISYM,
                        &ChemicalEvolutionTestNumberOfStars);
    ret += sscanf(line, "ChemicalEvolutionTestGasDensity = %"FSYM,
                        &ChemicalEvolutionTestGasDensity);
    ret += sscanf(line, "ChemicalEvolutionTestBackgroundGasDensity = %"FSYM,
                        &ChemicalEvolutionTestBackgroundGasDensity);
    ret += sscanf(line, "ChemicalEvolutionTestConcentration = %"FSYM,
                        &ChemicalEvolutionTestConcentration);
    ret += sscanf(line, "ChemicalEvolutionTestGasRadius = %"FSYM,
                        &ChemicalEvolutionTestGasRadius);

    ret += sscanf(line, "ChemicalEvolutionTestGasDistribution = %"ISYM,
                        &ChemicalEvolutionTestGasDistribution); // 0: uniform, 1: spherical isothermal

    ret += sscanf(line, "ChemicalEvolutionTestGasTemperature = %"FSYM,
                        &ChemicalEvolutionTestGasTemperature);
    ret += sscanf(line, "ChemicalEvolutionTestBackgroundGasTemperature = %"FSYM,
                        &ChemicalEvolutionTestBackgroundGasTemperature);

    ret += sscanf(line, "ChemicalEvolutionTestGasMetallicity = %"FSYM,
                        &ChemicalEvolutionTestGasMetallicity);
    ret += sscanf(line, "ChemicalEvolutionTestRefineAtStart = %"ISYM,
                        &ChemicalEvolutionTestRefineAtStart);

    ret += sscanf(line, "ChemicalEvolutionTestUseMetals = %"ISYM,
                        &ChemicalEvolutionTestUseMetals);

    ret += sscanf(line, "ChemicalEvolutionTestHFraction = %"FSYM,
                        &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "ChemicalEvolutionTestHIIFraction = %"FSYM,
                        &TestProblemData.HII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestHeIIFraction = %"FSYM,
                        &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestHeIIIFraction = %"FSYM,
                        &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestHMFraction = %"FSYM,
                        &TestProblemData.HM_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestH2IFraction = %"FSYM,
                        &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestH2IIFraction = %"FSYM,
                        &TestProblemData.H2II_Fraction);
    ret += sscanf(line, "ChemicalEvolutionTestDeuteriumToHydrogenRatio = %"FSYM,
                        &TestProblemData.DeuteriumToHydrogenRatio);

    ret += sscanf(line, "ChemicalEvolutionTestSpeciesFractions = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 0,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 1,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 2,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 3,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 4,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 5,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 6,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 7,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 8,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 9);

    ret += sscanf(line, "ChemicalEvolutionTestScaledSolarAbundances = %"ISYM,
                        &ChemicalEvolutionTestScaledSolarAbundances);

    ret += sscanf(line, "ChemicalEvolutionTestStarPosition = %"PSYM" %"PSYM" %"PSYM,
                       ChemicalEvolutionTestStarPosition, ChemicalEvolutionTestStarPosition+1, ChemicalEvolutionTestStarPosition+2);
    ret += sscanf(line, "ChemicalEvolutionTestStarVelocity = %"PSYM" %"PSYM" %"PSYM,
                     ChemicalEvolutionTestStarVelocity, ChemicalEvolutionTestStarVelocity+1,ChemicalEvolutionTestStarVelocity+2);
    ret += sscanf(line, "ChemicalEvolutionTestStarMass = %"FSYM,
                       &ChemicalEvolutionTestStarMass);
    ret += sscanf(line, "ChemicalEvolutionTestStarMetallicity = %"FSYM,
                       &ChemicalEvolutionTestStarMetallicity);
    ret += sscanf(line, "ChemicalEvolutionTestStarLifetime = %"FSYM,
                       &ChemicalEvolutionTestStarLifetime);


    if (ret == 0 && strstr(line, "=") && strstr(line, "ChemicalEvolutionTest")
                 && line[0] != '#' && !strstr(line, "ChemicalEvolutionTestStar")){
      fprintf(stderr, "warning: the following parameter line was not interpretd:\n%s\n",line);
    }
  }

  TestProblemData.MultiSpecies = MultiSpecies;
  TestProblemData.UseMetallicityField = ChemicalEvolutionTestUseMetals;
  /* check units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;

  if(GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
              &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  ChemicalEvolutionTestGasDensity     /= DensityUnits;
  ChemicalEvolutionTestBackgroundGasDensity /= DensityUnits;
  ChemicalEvolutionTestGasTemperature /= TemperatureUnits;
  ChemicalEvolutionTestBackgroundGasTemperature /= TemperatureUnits;

  if (ExternalGravity == 30){
    InitializeDoublePowerDarkMatter();

    // this currently does not work for some reason.... possibly
    // due to problems in computing acceleration from potential? idk...
    // externalgravity = 1 computes acceleration directly (analytically)...
    /// probably better anyway....

    if (DiskGravityDarkMatterDensity < 0){
        const float rho_crit = 9.33E-30; // cgs - make consistent in Grid_ChemicalEvolutionTest
        const float conc = ChemicalEvolutionTestConcentration; // temp
        DiskGravityDarkMatterDensity = 200.0 / 3.0 * rho_crit * conc*conc*conc / (log(1.0+conc)+conc/(1.0+conc)); //cgs
        const float M200 = DiskGravityDarkMatterMassInterior * SolarMass; // cgs
        DiskGravityDarkMatterCutoffR = POW((3.0*M200/(4.0*pi*200.0*rho_crit)),1.0/3.0) / Mpc_cm;
        DiskGravityDarkMatterR = DiskGravityDarkMatterCutoffR / conc;
    }// end if density

  } else if (ExternalGravity == 1){

    const float rho_crit = 9.33E-30; // cgs - make consistent in Grid_ChemicalEvolutionTest
    const float conc = ChemicalEvolutionTestConcentration;
    const float M200 = DiskGravityDarkMatterMassInterior * SolarMass; // cgs
    DiskGravityDarkMatterCutoffR = POW((3.0*M200/(4.0*pi*200.0*rho_crit)),1.0/3.0) / Mpc_cm;
    DiskGravityDarkMatterR = DiskGravityDarkMatterCutoffR / conc;
    DiskGravityDarkMatterDensity = 200.0 / 3.0 * rho_crit * conc*conc*conc / (log(1.0+conc)+conc/(1.0+conc)); //cgs

    HaloCentralDensity = DiskGravityDarkMatterDensity; // must be cgs
    HaloConcentration  = ChemicalEvolutionTestConcentration;
    HaloVirialRadius   = DiskGravityDarkMatterCutoffR * Mpc_cm; // cgs
  }

  // initialize star properties
  IndividualStarProperties_Initialize(MetaData);
  IndividualStarRadiationProperties_Initialize();
  InitializeStellarYields(MetaData.Time);



  /* set up grid */
  float BackgroundGasDensity     = ChemicalEvolutionTestGasDensity;
  float BackgroundGasTemperature = ChemicalEvolutionTestGasTemperature;

  if (ChemicalEvolutionTestGasDistribution == 1){
    BackgroundGasDensity = ChemicalEvolutionTestBackgroundGasDensity;
    BackgroundGasTemperature = ChemicalEvolutionTestBackgroundGasTemperature;
  } else if (ChemicalEvolutionTestGasDistribution > 1){
    ENZO_FAIL("ChemicalEvolutionTest: incorrect gas distribution type.");
  }

  float uniform_velocity[3]  = {0.0, 0.0, 0.0};
  float uniform_B_field[3]   = {0.0, 0.0, 0.0};
  float uniform_total_energy = BackgroundGasTemperature / ((Gamma-1.0)*0.6);

  if (TopGrid.GridData->InitializeUniformGrid(BackgroundGasDensity,
                                              uniform_total_energy,
                                              uniform_total_energy,
                                              uniform_velocity,
                                              uniform_B_field) == FAIL){
      ENZO_FAIL("Error in InitializeUniformGrid.");
  }

  if (TopGrid.GridData->ChemicalEvolutionTestInitializeGrid(ChemicalEvolutionTestGasDensity,
                                                            ChemicalEvolutionTestGasTemperature,
                                                            ChemicalEvolutionTestGasMetallicity,
                                                            TRUE) == FAIL){
    ENZO_FAIL("Error in ChemicalEvolutionInitialize[Sub]Grid.");
  } // end subgrid if

  /* add refinement levels */

  if (ChemicalEvolutionTestRefineAtStart){
    /* Initialize level array and fill */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

    for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++){
      LevelArray[level] = NULL;
    }

    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels tothe maximum depth or until no new levels are made */
    /* Restart each level once created */

    for (int level = 0; level < MaximumRefinementLevel; level++){
      printf("In level %"ISYM"\n", level);

      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL){
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }

      if (LevelArray[level+1] == NULL){
        break;
      }

      LevelHierarchyEntry *TempGrid = LevelArray[level+1];

      while (TempGrid != NULL) {
        TempGrid->GridData->InitializeUniformGrid(BackgroundGasDensity,
                                                  uniform_total_energy,
                                                  uniform_total_energy,
                                                  uniform_velocity,
                                                  uniform_B_field);

        TempGrid->GridData->ChemicalEvolutionTestInitializeGrid(ChemicalEvolutionTestGasDensity,
                                                                ChemicalEvolutionTestGasTemperature,
                                                                ChemicalEvolutionTestGasMetallicity,
                                                                FALSE);
        TempGrid = TempGrid->NextGridThisLevel;
      }
    } // end loop over levels

    /* Loop back from the bottom */
    for (int level = MaximumRefinementLevel; level > 0 ; level--){
      LevelHierarchyEntry *TempGrid = LevelArray[level];

      while (TempGrid != NULL){
        if (TempGrid->GridData->ProjectSolutionToParentGrid(
                                         *LevelArray[level-1]->GridData) == FAIL){
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
        TempGrid = TempGrid->NextGridThisLevel;
      }// end while
    } // end loop over levels
  } // end refine


  /* set up output field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  if(MetaData.TopGridRank > 1)
    DataLabel[count++] = Vel2Name;
  if(MetaData.TopGridRank > 2)
    DataLabel[count++] = Vel3Name;
  if(CRModel)
    DataLabel[count++] = CRName;

  if(WritePotential)
    DataLabel[count++] = GravPotentialName;

  /* handle the multispecies things */
  if (TestProblemData.MultiSpecies) {
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1){
      DataLabel[count++] = HMName;
      DataLabel[count++] = H2IName;
      DataLabel[count++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2){
      DataLabel[count++] = DIName;
      DataLabel[count++] = DIIName;
      DataLabel[count++] = HDIName;
    }
  }

  /* Metallicity and Metals */
  if (TestProblemData.UseMetallicityField){
    DataLabel[count++] = MetalName;

    if(MultiMetals ==2){

      for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        if(StellarYieldsAtomicNumbers[i] > 2){
          DataLabel[count++] = ChemicalSpeciesBaryonFieldLabel(StellarYieldsAtomicNumbers[i]);
        }
      }
    }

    if (IndividualStarTrackAGBMetalDensity){
      DataLabel[count++] = AGBMetalName;
    }

    if (IndividualStarPopIIIFormation){
      DataLabel[count++] = PopIIIMetalName;
      DataLabel[count++] = PopIIIPISNeMetalName;

      if (IndividualStarPopIIISeparateYields){
        for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
          if(StellarYieldsAtomicNumbers[i] > 2){
            DataLabel[count++] = ChemicalSpeciesBaryonFieldLabel(StellarYieldsAtomicNumbers[i],2);
          }
        }
      }
    }

    if (IndividualStarTrackWindDensity){
      DataLabel[count++] = WindMetalName;
      DataLabel[count++] = WindMetalName2;
    }

    if (IndividualStarTrackSNMetalDensity){
      DataLabel[count++] = SNIaMetalName;
      if (IndividualStarSNIaModel == 2){
        DataLabel[count++] = ExtraMetalName0;
        DataLabel[count++] = ExtraMetalName1;
        DataLabel[count++] = ExtraMetalName2;
      }
      DataLabel[count++] = SNIIMetalName;
    }

    if (IndividualStarRProcessModel){
      DataLabel[count++] = RProcMetalName;
    }
  }

 // fill in remaining slots
 for(int j=0; j < count; j++) DataUnits[j] = NULL;

 /* Write Parameters to parameter output file */
 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "ChemicalEvolutionTestGasDensity = %"GSYM"\n", ChemicalEvolutionTestGasDensity);
   fprintf(Outfptr, "ChemicalEvolutionTestGasTemperature = %"GSYM"\n", ChemicalEvolutionTestGasTemperature);
   fprintf(Outfptr, "ChemicalEvolutionTestGasMetallicity = %"FSYM"\n", ChemicalEvolutionTestGasMetallicity);

   fprintf(Outfptr, "ChemicalEvolutionTestBackgroundGasDensity = %"GSYM"\n", ChemicalEvolutionTestBackgroundGasDensity);
   fprintf(Outfptr, "ChemicalEvolutionTestBackgroundGasTemperature = %"GSYM"\n", ChemicalEvolutionTestBackgroundGasTemperature);

   fprintf(Outfptr, "ChemicalEvolutionTestRefineAtStart = %"ISYM"\n", ChemicalEvolutionTestRefineAtStart);

   fprintf(Outfptr, "CHemicalEvolutionTestUseMetals = %"ISYM"\n", ChemicalEvolutionTestUseMetals);

   fprintf(fptr, "ChemicalEvolutionTestStarPosition = ");
   WriteListOfFloats(fptr, MetaData.TopGridRank, ChemicalEvolutionTestStarPosition);
   fprintf(fptr, "ChemicalEvolutionTestStarVelocity = ");
   WriteListOfFloats(fptr, MetaData.TopGridRank, ChemicalEvolutionTestStarVelocity);
   fprintf(fptr, "ChemicalEvolutionTestStarMass = %"FSYM"\n", ChemicalEvolutionTestStarMass);
   fprintf(fptr, "ChemicalEvolutionTestStarMetallicity = %"FSYM"\n", ChemicalEvolutionTestStarMetallicity);
   fprintf(fptr, "ChemicalEvolutionTestStarLifetime = %"FSYM"\n", ChemicalEvolutionTestStarLifetime);
   fprintf(fptr, "ChemicalEvolutionTestScaledSolarAbundances = %"ISYM"\n", ChemicalEvolutionTestScaledSolarAbundances);


   fprintf(Outfptr, "ChemicalEvolutionTestHydryogenFractionByMass = %"FSYM"\n", TestProblemData.HydrogenFractionByMass);
   fprintf(Outfptr, "ChemicalEvolutionTestHIIFraction = %"FSYM"\n", TestProblemData.HII_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestHeIIFraction = %"FSYM"\n", TestProblemData.HeII_Fraction);
   fprintf(Outfptr, "ChemicalEvolutionTestHeIIIFraction = %"FSYM"\n", TestProblemData.HeIII_Fraction);
   fprintf(Outfptr, "ChemicalEvolutionTestHMFraction = %"FSYM"\n", TestProblemData.HM_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestH2IFraction = %"FSYM"\n", TestProblemData.H2I_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestH2IIFraction = %"FSYM"\n", TestProblemData.H2II_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio );

   fprintf(Outfptr, "ChemicalEvolutionTestCIFraction = %"FSYM"\n", TestProblemData.CI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestOIFraction = %"FSYM"\n", TestProblemData.OI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestNIFraction = %"FSYM"\n", TestProblemData.NI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestMgIFraction = %"FSYM"\n", TestProblemData.MgI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestSiIFraction = %"FSYM"\n", TestProblemData.SiI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestFeIFraction = %"FSYM"\n", TestProblemData.FeI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestLaIFraction = %"FSYM"\n", TestProblemData.LaI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestBaIFraction = %"FSYM"\n", TestProblemData.BaI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestYIFraction = %"FSYM"\n", TestProblemData.YI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionTestEuIFraction = %"FSYM"\n", TestProblemData.EuI_Fraction );



   fprintf(Outfptr, "ChemicalEvolutionTestSpeciesFractions   = ");
   WriteListOfFloats(Outfptr, MAX_STELLAR_YIELDS, TestProblemData.ChemicalTracerSpecies_Fractions);

 }



}
