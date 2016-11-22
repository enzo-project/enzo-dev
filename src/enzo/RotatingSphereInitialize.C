/***********************************************************************
/
/   Rotating Sphere Test Problem
/
/   written by: Greg Meece, Brian O'Shea
/   date:          July 2013
/   modified1:   March 2014 to clean up the code.
/
/   PURPOSE: Sets up a rotating sphere of gas in an NFW potential.
/                Originally written to simulate Pop III Star Formation.
/                For details, see Meece (2014) in ApJ.
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
 
int RotatingSphereInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
          TopGridData &MetaData)
{
   if(debug){
      printf("Entering RotatingSphereInitialize\n");
      fflush(stdout);
      }

   // Make sure that we are working in 3D
   if (MetaData.TopGridRank != 3) {
      ENZO_VFAIL("Cannot do RotatingSphere in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
      }

   // Field Names
   char *DensName = "Density";
   char *TEName    = "TotalEnergy";
   char *GEName    = "GasEnergy";
   char *Vel1Name = "x-velocity";
   char *Vel2Name = "y-velocity";
   char *Vel3Name = "z-velocity";
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
 
   FLOAT RotatingSphereSubgridLeft[MAX_DIMENSION], RotatingSphereSubgridRight[MAX_DIMENSION];
   FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

   float RotatingSphereNFWMass;
   float RotatingSphereNFWConcentration;
   float RotatingSphereCoreRadius;
   float RotatingSphereCentralDensity;
   float RotatingSphereCoreDensityExponent;
   float RotatingSphereOuterDensityExponent;
   float RotatingSphereExternalTemperature;
   float RotatingSphereSpinParameter;
   float RotatingSphereAngularMomentumExponent;
   int RotatingSphereUseTurbulence;
   float RotatingSphereTurbulenceRMS;
   float RotatingSphereRedshift;

   // Local variable declarations
   char line[MAX_LINE_LENGTH];
   int   i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
      SubgridDims[MAX_DIMENSION];

   // Initialize parameters to default values
   RotatingSphereNFWMass = 1.0e7;
   RotatingSphereNFWConcentration = 2;
   RotatingSphereCoreRadius = 16.0;
   RotatingSphereCentralDensity = 1.0;
   RotatingSphereCoreDensityExponent = 0.1;
   RotatingSphereOuterDensityExponent = 2.5;
   RotatingSphereExternalTemperature = 200.0;
   RotatingSphereSpinParameter = 0.05;
   RotatingSphereAngularMomentumExponent = 0.9;
   RotatingSphereUseTurbulence = 1;
   RotatingSphereTurbulenceRMS = 0.01;
   RotatingSphereRedshift = 20.0;

   // Set no subgrids by default. 
   RotatingSphereSubgridLeft[0] = RotatingSphereSubgridLeft[1] = 
      RotatingSphereSubgridLeft[2] = 0.0;      // start of subgrid(s)

   RotatingSphereSubgridRight[0] = RotatingSphereSubgridRight[1] = 
      RotatingSphereSubgridRight[2] = 0.0;      // end of subgrid(s)

   // Read problem specific parameters. 
   while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
      ret = 0;

      ret += sscanf(line, "RotatingSphereNFWMass = %"FSYM, &RotatingSphereNFWMass);
      ret += sscanf(line, "RotatingSphereNFWConcentration = %"FSYM, &RotatingSphereNFWConcentration);
      ret += sscanf(line, "RotatingSphereCoreRadius = %"FSYM, &RotatingSphereCoreRadius);
      ret += sscanf(line, "RotatingSphereCentralDensity = %"FSYM, &RotatingSphereCentralDensity);
      ret += sscanf(line, "RotatingSphereCoreDensityExponent = %"FSYM, &RotatingSphereCoreDensityExponent);
      ret += sscanf(line, "RotatingSphereOuterDensityExponent = %"FSYM, &RotatingSphereOuterDensityExponent);
      ret += sscanf(line, "RotatingSphereExternalTemperature = %"FSYM, &RotatingSphereExternalTemperature);
      ret += sscanf(line, "RotatingSphereSpinParameter = %"FSYM, &RotatingSphereSpinParameter);
      ret += sscanf(line, "RotatingSphereAngularMomentumExponent = %"FSYM, &RotatingSphereAngularMomentumExponent);
      ret += sscanf(line, "RotatingSphereUseTurbulence = %"ISYM, &RotatingSphereUseTurbulence);
      ret += sscanf(line, "RotatingSphereTurbulenceRMS = %"FSYM, &RotatingSphereTurbulenceRMS);
      ret += sscanf(line, "RotatingSphereRedshift = %"FSYM, &RotatingSphereRedshift);

      ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);
      ret += sscanf(line, "TestProblemDeuteriumToHydrogenRatio = %"FSYM, &TestProblemData.DeuteriumToHydrogenRatio);

      ret += sscanf(line, "TestProblemInitialHIFraction = %"FSYM, &TestProblemData.HI_Fraction);
      ret += sscanf(line, "TestProblemInitialHIIFraction = %"FSYM, &TestProblemData.HII_Fraction);
      ret += sscanf(line, "TestProblemInitialHeIFraction = %"FSYM, &TestProblemData.HeI_Fraction);
      ret += sscanf(line, "TestProblemInitialHeIIFraction = %"FSYM, &TestProblemData.HeII_Fraction);
      ret += sscanf(line, "TestProblemInitialHeIIIFraction = %"FSYM, &TestProblemData.HeIII_Fraction);
      ret += sscanf(line, "TestProblemInitialHMFraction = %"FSYM, &TestProblemData.HM_Fraction);
      ret += sscanf(line, "TestProblemInitialH2IFraction = %"FSYM, &TestProblemData.H2I_Fraction);
      ret += sscanf(line, "TestProblemInitialH2IIFraction = %"FSYM, &TestProblemData.H2II_Fraction);
      ret += sscanf(line, "TestProblemInitialDIFraction = %"FSYM, &TestProblemData.DI_Fraction);
      ret += sscanf(line, "TestProblemInitialDIIFraction = %"FSYM, &TestProblemData.DII_Fraction);
      ret += sscanf(line, "TestProblemInitialHDIFraction = %"FSYM, &TestProblemData.HDI_Fraction);

      ret += sscanf(line, "TestProblemUseMetallicityField = %"ISYM, &TestProblemData.UseMetallicityField);
      ret += sscanf(line, "TestProblemInitialMetallicityFraction = %"FSYM, &TestProblemData.MetallicityField_Fraction);


      // Issue a warning if the line is suspicious 
      if (ret == 0 && strstr(line, "=") && (strstr(line, "RotatingSphere") || strstr(line, "TestProblem")) &&
   line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
         fprintf(stderr,
                 "*** warning: the following parameter line was not interpreted:\n%s\n",
                 line);
      } // end input from parameter file

   TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)
 
   // Initialize a uniform grid
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
 
   /* Create as many subgrids as there are refinement levels
       needed to resolve the initial explosion region upon the start-up. */
 
//   HierarchyEntry ** Subgrid;
//   if (MaximumRefinementLevel > 0)
//      Subgrid    = new HierarchyEntry*[MaximumRefinementLevel];
// 
//   /* Create new HierarchyEntries. */
// 
//   int lev;
//   for (lev = 0; lev < MaximumRefinementLevel; lev++)
//      Subgrid[lev] = new HierarchyEntry;
// 
//   for (lev = 0; lev < MaximumRefinementLevel; lev++) {
// 
//      for (dim = 0; dim < MetaData.TopGridRank; dim++)
//         NumberOfSubgridZones[dim] =
//   nint((RotatingSphereSubgridRight[dim] - RotatingSphereSubgridLeft[dim])/
//          ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
//            float(MetaData.TopGridDims[dim])))
//            *int(POW(RefineBy, lev + 1));
// 
//      if (debug)
//         printf("RotatingSphere:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
//          NumberOfSubgridZones[0]);
// 
//      if (NumberOfSubgridZones[0] > 0) {
// 
//         /* fill them out */
// 
//         if (lev == 0)
//   TopGrid.NextGridNextLevel   = Subgrid[0];
//         Subgrid[lev]->NextGridThisLevel = NULL;
//         if (lev == MaximumRefinementLevel-1)
//   Subgrid[lev]->NextGridNextLevel = NULL;
//         else
//   Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
//         if (lev == 0)
//   Subgrid[lev]->ParentGrid            = &TopGrid;
//         else
//   Subgrid[lev]->ParentGrid            = Subgrid[lev-1];
// 
//         /* compute the dimensions and left/right edges for the subgrid */
// 
//         for (dim = 0; dim < MetaData.TopGridRank; dim++) {
//   SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
//   LeftEdge[dim]      = RotatingSphereSubgridLeft[dim];
//   RightEdge[dim]    = RotatingSphereSubgridRight[dim];
//         }
// 
//         /* create a new subgrid and initialize it */
// 
//         Subgrid[lev]->GridData = new grid;
//         Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
//         Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
//                   LeftEdge, RightEdge, 0);
//         if (Subgrid[lev]->GridData->InitializeUniformGrid(uniform_density,
//                      uniform_total_energy,
//                      uniform_total_energy,
//                      uniform_velocity,
//                      uniform_B_field) == FAIL) {
//      ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
//         }
// 
//         /* set up the initial explosion area on the finest resolution subgrid */
// 
//         if (lev == MaximumRefinementLevel - 1)
//   if (Subgrid[lev]->GridData->RotatingSphereInitializeGrid(RotatingSphereNFWMass,
//                               RotatingSphereNFWConcentration,
//                               RotatingSphereCoreRadius,
//                               RotatingSphereCentralDensity,
//                               RotatingSphereCoreDensityExponent,
//                               RotatingSphereOuterDensityExponent,
//                               RotatingSphereExternalTemperature,
//                               RotatingSphereSpinParameter,
//                               RotatingSphereAngularMomentumExponent,
//                               RotatingSphereUseTurbulence,
//                               RotatingSphereTurbulenceRMS,
//                               RotatingSphereRedshift)
//         == FAIL) {
//            ENZO_FAIL("Error in RotatingSphereInitialize[Sub]Grid.");
//   }
//
//      }
//      else{
//         printf("RotatingSphere: single grid start-up.\n");
//      }
//   }
//
// 
//   /* set up subgrids from level 1 to max refinement level -1 */
// 
//   for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
//      if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
//                      *(Subgrid[lev-1]->GridData))
//   == FAIL) {
//                  ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
//      }
// 
//   /* set up the root grid */
// 
//   if (MaximumRefinementLevel > 0) {
//      if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
//   == FAIL) {
//                  ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
//      }
//   }
//   else
      if (TopGrid.GridData->RotatingSphereInitializeGrid(RotatingSphereNFWMass,
                               RotatingSphereNFWConcentration,
                               RotatingSphereCoreRadius,
                               RotatingSphereCentralDensity,
                               RotatingSphereCoreDensityExponent,
                               RotatingSphereOuterDensityExponent,
                               RotatingSphereExternalTemperature,
                               RotatingSphereSpinParameter,
                               RotatingSphereAngularMomentumExponent,
                               RotatingSphereUseTurbulence,
                               RotatingSphereTurbulenceRMS,
                               RotatingSphereRedshift) ==FAIL) {
                  ENZO_FAIL("Error in RotatingSphereInitializeGrid.");
      }

   // If desired, refine the grid during initialization
   int RotatingSphereRefineAtStart = 1;

   if (RotatingSphereRefineAtStart) {
      /* Declare, initialize and fill out the LevelArray. */
 
      LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

      for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
         LevelArray[level] = NULL;

      AddLevel(LevelArray, &TopGrid, 0);
 
      /* Add levels to the maximum depth or until no new levels are created,
         and re-initialize the level after it is created. */
 
      for (int level = 0; level < MaximumRefinementLevel; level++) {
         printf("In level %"ISYM"\n", level);

         if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
            fprintf(stderr, "Error in RebuildHierarchy.\n");
            return FAIL;
            }

         if (LevelArray[level+1] == NULL)
            break;
 
         printf("Going to create a new grid!\n");
         fflush(stdout);
 
         LevelHierarchyEntry *Temp = LevelArray[level+1];
         while (Temp != NULL) {
            Temp->GridData->RotatingSphereInitializeGrid(RotatingSphereNFWMass,
                                                       RotatingSphereNFWConcentration,
                                                       RotatingSphereCoreRadius,
                                                       RotatingSphereCentralDensity,
                                                       RotatingSphereCoreDensityExponent,
                                                       RotatingSphereOuterDensityExponent,
                                                       RotatingSphereExternalTemperature,
                                                       RotatingSphereSpinParameter,
                                                       RotatingSphereAngularMomentumExponent,
                                                       RotatingSphereUseTurbulence,
                                                       RotatingSphereTurbulenceRMS,
                                                       RotatingSphereRedshift);
            Temp = Temp->NextGridThisLevel;
            }
         } // end: loop over levels
 
 
      /* Loop back from the bottom, restoring the consistency among levels. */
 
      for (int level = MaximumRefinementLevel; level > 0; level--) {
         LevelHierarchyEntry *Temp = LevelArray[level];

         while (Temp != NULL) {
            if (Temp->GridData->ProjectSolutionToParentGrid(
         			    *LevelArray[level-1]->GridData) == FAIL) {
               fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
               return FAIL;
               }

            Temp = Temp->NextGridThisLevel;
            }
         }
      } // end: if (RotatingSphereRefineAtStart)

 
   /* set up field names and units -- NOTE: these absolutely MUST be in 
       the same order that they are in Grid_InitializeUniformGrids.C, or 
       else you'll find out that data gets written into incorrectly-named
       fields.   Just FYI. */

   i = 0;
   DataLabel[i++] = DensName;
   DataLabel[i++] = TEName;
   if(DualEnergyFormalism)
      DataLabel[i++] = GEName;
   DataLabel[i++] = Vel1Name;

   if(MetaData.TopGridRank > 1)
      DataLabel[i++] = Vel2Name;

   if(MetaData.TopGridRank > 2)
      DataLabel[i++] = Vel3Name;

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
 
   if(debug){
      printf("Exiting RotatingSphereInitialize\n");
      fflush(stdout);
      }
 
   return SUCCESS;
}
