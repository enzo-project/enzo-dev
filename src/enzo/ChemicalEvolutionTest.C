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

  /* Names for chemical evolution element abundances */
  char *CIName  =  "CI_Density";
  char *OIName  =  "OI_Density";
  char *NIName  =  "NI_Density";
  char *MgIName = "MgI_Density";
  char *SiIName = "SiI_Density";
  char *FeIName = "FeI_Density";
  char *YIName  =  "YI_Density";
  char *LaIName = "LaI_Density";
  char *BaIName = "BaI_Density";
  char *EuIName = "EuI_Density";
   
  /* declarations */
  char line[MAX_LINE_LENGTH];
  int  dim, ret, level, disk, i; // might not need disk

  /* make sure we are in 3D */
  
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do ChemicalEvolutionTest in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */
  float GasDensity     = 1.673E-24,
        GasTemperature = 1.0E4,
        StarPosX       = 0.0,
        StarPosY       = 0.0,
        StarPosZ       = 0.0,
        StarMass       = 1.61109E34 ; // 8.1 solar masses
         // need to change defaults before ready 1/21/16
 

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
   
    ret = 0;

    ret += sscanf(line, "ChemicalEvolutionGasDensity = %"FSYM, &GasDensity);
    ret += sscanf(line, "ChemicalEvolutionGasTemperature = %"FSYM, &GasTemperature);
    ret += sscanf(line, "ChemicalEvolutionStarPosX = %"FSYM, &StarPosX);
    ret += sscanf(line, "ChemicalEvolutionStarPosY = %"FSYM, &StarPosY);
    ret += sscanf(line, "ChemicalEvolutionStarPosZ = %"FSYM, &StarPosZ);
    ret += sscanf(line, "ChemicalEvolutionStarMass = %"FSYM, &StarMass);

    ret += sscanf(line, "ChemicalEvolutionInitialHIFraction = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialHIIFraction = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialHeIFraction = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialHeIIFraction = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialHeIIIFraction = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialHMFraction = %"FSYM, &TestProblemData.HM_Fraction);

    ret += sscanf(line, "ChemicalEvolutionInitialDIFraction  = %"FSYM, &TestProblemData.DI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionProblemInitialDIIFraction  = %"FSYM, &TestProblemData.DII_Fraction);
    ret += sscanf(line, "ChemicalEvolutionProblemInitialHDIFraction  = %"FSYM, &TestProblemData.HDI_Fraction);

    ret += sscanf(line, "ChemicalEvolutionInitialCIFraction = %"FSYM, &TestProblemData.CI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialNIFraction = %"FSYM, &TestProblemData.NI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialOIFraction = %"FSYM, &TestProblemData.OI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialMgIFraction = %"FSYM, &TestProblemData.MgI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialSiIFraction = %"FSYM, &TestProblemData.SiI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialFeIFraction = %"FSYM, &TestProblemData.FeI_Fracrtion);
    ret += sscanf(line, "ChemicalEvolutionInitialYIFraction = %"FSYM, &TestProblemData.YI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialLaIFraction = %"FSYM, &TestProblemData.LaI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialBaIFraction = %"FSYM, &TestProblemData.BaI_Fraction);
    ret += sscanf(line, "ChemicalEvolutionInitialEuIFraction = %"FSYM, &TestProblemData.EuI_Fraction);
  }

  /* check units */
  float DensityUnits, LengthUnits, TemperatureUnites, TimeUnits, VelocityUnits, MassUnits;

  if(GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
              &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  GasDensity = GasDensity / DensityUnits;
  GasTemperature = GasTemperature / TemperatureUnits;
  StarMass = StarMass / MassUnits;
  StarPosX = StarPosX / LengthUnits; 
  StarPosY = StarPosY / LengthUnits;
  StarPosZ = StarPosZ / LengthUnits;

  /* set up grid */

  if (TopGrid.GridData->ChemicalEvolutionInitializeGrid(GasDensity, GasTemperature,
                                                        StarPosX, StarPoxY, StarPosZ,
                                                        StarMass) == FAIL){
    ENZO_FAIL("Error in ChemicalEvolutionInitialize[Sub]Grid.");
  } // end subgrid if

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
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
   
    if (TestProblemData.MultiSpecies > 1){
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    } 
    if (TestProblemData.MultiSpecies > 2){
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }

  /* Metallicity and Metals */
  if (TestProblemData.UseMetallicityField){
    DataLabel[i++] = MetalName;

    int MM = TestProblemData.MultiMetals; // for convenience

    if( MM == 1){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
    if( (MM==2) || (MM==5) || (MM==6) || (MM==9)){
      DataLabel[i++] = CIName;
      DataLabel[i++] = OIName;
      DataLabel[i++] = NIName;
      DataLabel[i++] = MgIName;
      DataLabel[i++] = SiIName;
      DataLabel[i++] = FeIName;
    }
    if( (MM==3) || (MM==5) || (MM==7) || (MM==9)){
      DataLabel[i++] = YIName;
      DataLabel[i++] = BaIName;
      DataLabel[i++] = LaIName;
    }


 // fill in remaining slots
 for(j=0; j < i; j++) DataUnits[j] = NULL;

 /* Write Parameters to parameter output file */
 if (MyProcessorNumber == ROOT_PROCESSOR) {
   fprintf(Outfptr, "ChemicalEvolutionInitialHIFraction = %"FSYM"\n", TestProblemData.HI_Fraction);
   fprintf(Outfptr, "ChemicalEvolutionInitialHIIFraction = %"FSYM"\n", TestProblemData.HII_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialHeIFraction = %"FSYM"\n", TestProblemData.HeI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialHeIIFraction = %"FSYM"\n", TestProblemData.HeII_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialHMFraction = %"FSYM"\n", TestProblemData.HM_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialDIFraction = %"FSYM"\n", TestProblemData.DI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialDIIFraction = %"FSYM"\n", TestProblemData.DII_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialHDIFraction = %"FSYM"\n", TestProblemData.HDI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialCIFraction = %"FSYM"\n", TestProblemData.CI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialOIFraction = %"FSYM"\n", TestProblemData.OI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialNIFraction = %"FSYM"\n", TestProblemData.NI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialMgIFraction = %"FSYM"\n", TestProblemData.MgI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialSiIFraction = %"FSYM"\n", TestProblemData.SiI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialFeIFraction = %"FSYM"\n", TestProblemData.FeI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialLaIFraction = %"FSYM"\n", TestProblemData.LaI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialBaIFraction = %"FSYM"\n", TestProblemData.BaI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialYIFraction = %"FSYM"\n", TestProblemData.YI_Fraction );
   fprintf(Outfptr, "ChemicalEvolutionInitialEuIFraction = %"FSYM"\n", TestProblemData.EuI_Fraction );
 
 } 

 
}


