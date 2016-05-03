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
#include "TopGridData.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level);

char* ChemicalSpeciesBaryonFieldLabel(const int &atomic_number);

void WriteListOfFloats(FILE *fptr, int N, float floats[]);

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
  char *MetalName   = "Metal_Density";
  char *MetalIaName = "MetalSNIa_Density";


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

  char *PeHeatingRateName = "Pe_heating_rate";


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
  float GasDensity     = 1.673E-24,
        GasTemperature = 1.0E4,
        GasMetallicity = tiny_number;

  int   ChemicalEvolutionTestRefineAtStart  = 1,
        ChemicalEvolutionTestUseMetals      = 1;

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

    ret += sscanf(line, "ChemicalEvolutionTestGasDensity = %"FSYM, 
                        &GasDensity);
    ret += sscanf(line, "ChemicalEvolutionTestGasTemperature = %"FSYM,
                        &GasTemperature);
    ret += sscanf(line, "ChemicalEvolutionTestGasMetallicity = %"FSYM,
                        &GasMetallicity);
    ret += sscanf(line, "ChemicalEvolutionTestRefineAtStart = %"ISYM,
                        &ChemicalEvolutionTestRefineAtStart);

    ret += sscanf(line, "ChemicalEvolutionTestMultiMetals = %"ISYM,
                        &TestProblemData.MultiMetals);
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

  GasDensity = GasDensity / DensityUnits;
  GasTemperature = GasTemperature / TemperatureUnits;

  /* set up grid */

  float uniform_velocity[3]  = {0.0, 0.0, 0.0};
  float uniform_B_field[3]   = {0.0, 0.0, 0.0};
  float uniform_total_energy = GasTemperature / ((Gamma-1.0)*0.6);

  if (TopGrid.GridData->InitializeUniformGrid(GasDensity,
                                              uniform_total_energy,
                                              uniform_total_energy,
                                              uniform_velocity,
                                              uniform_B_field) == FAIL){
      ENZO_FAIL("Error in InitializeUniformGrid.");
  }


  if (TopGrid.GridData->ChemicalEvolutionTestInitializeGrid(GasDensity,
                                                        GasTemperature,
                                                        GasMetallicity) == FAIL){
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
        TempGrid->GridData->InitializeUniformGrid(GasDensity,
                                                  uniform_total_energy,
                                                  uniform_total_energy,
                                                  uniform_velocity,
                                                  uniform_B_field);

        TempGrid->GridData->ChemicalEvolutionTestInitializeGrid(GasDensity,
                                                                GasTemperature,
                                                                GasMetallicity);
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

    if(TestProblemData.MultiMetals ==2){

      for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        if(StellarYieldsAtomicNumbers[i] > 2){
          DataLabel[count++] = ChemicalSpeciesBaryonFieldLabel(StellarYieldsAtomicNumbers[i]);
        }
      }
    }
  }

 if(STARMAKE_METHOD(INDIVIDUAL_STAR) && IndividualStarFUVHeating){
  DataLabel[count++] = PeHeatingRateName;
 }

 // fill in remaining slots
 for(int j=0; j < count; j++) DataUnits[j] = NULL;

 /* Write Parameters to parameter output file */
 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "ChemicalEvolutionTestGasDensity = %"GSYM"\n", GasDensity);
   fprintf(Outfptr, "ChemicalEvolutionTestGasTemperature = %"GSYM"\n", GasTemperature);
   fprintf(Outfptr, "ChemicalEvolutionTestGasMetallicity = %"FSYM"\n", GasMetallicity);
   fprintf(Outfptr, "ChemicalEvolutionTestRefineAtStart = %"ISYM"\n", ChemicalEvolutionTestRefineAtStart);

   fprintf(Outfptr, "ChemicalEvolutionTestMultiMetals = %"ISYM"\n", TestProblemData.MultiMetals);
   fprintf(Outfptr, "CHemicalEvolutionTestUseMetals = %"ISYM"\n", ChemicalEvolutionTestUseMetals);

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



