/***********************************************************************
/
/   Thermal Instability Test Problem
/
/   written by: Cameron Hummels, Iryna Butsky
/   date:       March 2018
/
/   PURPOSE: Investigate how thermal instability operates when different
/            physics are present
/
/   RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
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
#include "TopGridData.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
 
int ThermalInstabilityInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
          TopGridData &MetaData)
{
   if(debug){
      printf("Entering ThermalInstabilityInitialize\n");
      fflush(stdout);
      }

   // Make sure that we are working in 3D
   if (MetaData.TopGridRank != 3) {
      ENZO_VFAIL("Cannot do ThermalInstability in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
      }

   // Field Names
   char *DensName = "Density";
   char *TEName   = "TotalEnergy";
   char *GEName   = "GasEnergy";
   char *Vel1Name = "x-velocity";
   char *Vel2Name = "y-velocity";
   char *Vel3Name = "z-velocity";
   char *BxName   = "Bx";
   char *ByName = "By";
   char *BzName = "Bz";
   char *PhiName = "Phi";
   char *CRName = "CREnergyDensity";
   char *ElectronName = "Electron_Density";
   char *HIName = "HI_Density";
   char *HIIName = "HII_Density";
   char *HeIName = "HeI_Density";
   char *HeIIName = "HeII_Density";
   char *HeIIIName = "HeIII_Density";
   char *HMName = "HM_Density";
   char *H2IName = "H2I_Density";
   char *H2IIName = "H2II_Density";
   char *DIName = "DI_Density";
   char *DIIName = "DII_Density";
   char *HDIName = "HDI_Density";
   char *MetalName = "Metal_Density";

   /* parameter declarations */
 
   FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

   // Local variable declarations
   char line[MAX_LINE_LENGTH];
   int   dim, ret, level;

   // Initialize parameters to default values
   float TIMeanDensity = 1;
   float TIDensityPerturbationAmplitude = 0.01;
   float TIMeanTemperature = 1000000;

   // Read problem specific parameters. 
   while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
      ret = 0;

      ret += sscanf(line, "TIMeanDensity = %"FSYM, &TIMeanDensity);
      ret += sscanf(line, "TIDensityPerturbationAmplitude = %"FSYM, &TIDensityPerturbationAmplitude);
      ret += sscanf(line, "TIMeanTemperature = %"FSYM, &TIMeanTemperature);
      ret += sscanf(line, "TestProblemUseMetallicityField = %"ISYM, &TestProblemData.UseMetallicityField);

      // Issue a warning if the line is suspicious 
      if (ret == 0 && strstr(line, "=") && strstr(line, "TI") && line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
         fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    } // end input from parameter file

    TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

    /* set the periodic boundaries */

    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
        MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
        MetaData.RightFaceBoundaryCondition[dim] = periodic;
    }

    /* Initialize a uniform grid first */

   float uniform_velocity[3] = {0.0, 0.0, 0.0};
   float uniform_density = 1.0;
   float uniform_total_energy = 1.0;
   float uniform_B_field[3] = {0.0, 0.0, 0.0};

   if (TopGrid.GridData->InitializeUniformGrid(uniform_density,
                        uniform_total_energy,
                        uniform_total_energy,
                        uniform_velocity,
                        uniform_B_field) == FAIL) {
                        ENZO_FAIL("Error in InitializeUniformGrid.");
                        }
    
    /* set up grid fields from turbulent ICs*/
    if (TopGrid.GridData->ThermalInstabilityInitializeGrid(TIMeanDensity,
                                                           TIDensityPerturbationAmplitude,
                                                           TIMeanTemperature) == FAIL) {
        ENZO_FAIL("Error in ThermalInstabilityInitializeGrid.");
    }

    /* set up field names and units */

    int i = 0;
    int j;
    DataLabel[i++] = DensName;
    DataLabel[i++] = TEName;
    if (DualEnergyFormalism)
        DataLabel[i++] = GEName;
    DataLabel[i++] = Vel1Name;
    DataLabel[i++] = Vel2Name;
    DataLabel[i++] = Vel3Name;
    if (HydroMethod == MHD_RK) {
      DataLabel[i++] =  BxName;
      DataLabel[i++] =  ByName;
      DataLabel[i++] =  BzName;
      DataLabel[i++] =  PhiName;
    }
    if (CRModel)
      DataLabel[i++] = CRName;

    if (TestProblemData.MultiSpecies) {
      DataLabel[i++] = ElectronName;
      DataLabel[i++] = HIName;
      DataLabel[i++] = HIIName;
      DataLabel[i++] = HeIName;
      DataLabel[i++] = HeIIName;
      DataLabel[i++] = HeIIIName;

      if (TestProblemData.MultiSpecies > 1) {
         DataLabel[i++] = HMName;
         DataLabel[i++] = H2IName;
         DataLabel[i++] = H2IIName;
       }

      if (TestProblemData.MultiSpecies > 2) {
         DataLabel[i++] = DIName;
         DataLabel[i++] = DIIName;
         DataLabel[i++] = HDIName;
       }
    }
    if (TestProblemData.UseMetallicityField)
        DataLabel[i++] = MetalName;

   for(j=0; j < i; j++)
      DataUnits[j] = NULL;

    /* Write parameters to parameter output file */

    if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(Outfptr, "TIMeanDensity  = %"FSYM"\n", TIMeanDensity);
        fprintf(Outfptr, "TIDensityPerturbationAmplitude = %"FSYM"\n", TIDensityPerturbationAmplitude);
        fprintf(Outfptr, "TIMeanTemperature  = %"FSYM"\n", TIMeanTemperature);
    }

    return SUCCESS;

}
