/***********************************************************************
/
/  INITIALIZE A KELVIN-HELMHOLTZ INSTABILITY SIMULATION
/
/  written by: Greg Bryan
/  date:       February, 1995
/  modified1:  Alexei Kritsuk, December 2004.
/  modified2:  Gregg Dobrowalski, Feb 2005.
/  modified3:  Alexei Kritsuk, April 2005. added more parameters.
/  modified4:  Cameron Hummels, December 2013. Added ramp.
/
/  PURPOSE:
/    Periodic boundary conditions
/    at all boundaries.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

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
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
                     int level);
int KHInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                 TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";

  /* parameter declarations */

  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, level;

  /* set default parameters */

  int RefineAtStart             = TRUE; 
  float KHInnerPressure         = 2.5;
  float KHOuterPressure         = 2.5;
  float KHVelocityJump          = 1.0;
  float KHPerturbationAmplitude = 0.1;
  float KHInnerDensity          = 2.0;
  float KHOuterDensity          = 1.0;
  float KHBulkVelocity          = 0.0;
  int   KHRamp                  = 1;    // Convergent ICs with Ramp
  float KHRampWidth             = 0.05;
  int   KHRandomSeed            = 123456789;

  float KHInnerInternalEnergy, KHOuterInternalEnergy;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "KHInnerDensity  = %"FSYM, &KHInnerDensity);
    ret += sscanf(line, "KHOuterDensity  = %"FSYM, &KHOuterDensity);
    ret += sscanf(line, "KHInnerPressure = %"FSYM, &KHInnerPressure);
    ret += sscanf(line, "KHOuterPressure = %"FSYM, &KHOuterPressure);
    ret += sscanf(line, "KHVelocityJump  = %"FSYM, &KHVelocityJump);
    ret += sscanf(line, "KHPerturbationAmplitude = %"FSYM, 
                  &KHPerturbationAmplitude);
    ret += sscanf(line, "KHBulkVelocity  = %"FSYM, &KHBulkVelocity);
    ret += sscanf(line, "KHRamp = %"ISYM, &KHRamp);
    ret += sscanf(line, "KHRampWidth     = %"FSYM, &KHRampWidth);
    ret += sscanf(line, "KHRandomSeed    = %"ISYM, &KHRandomSeed);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "KH") && 
        line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
            "warning: the following parameter line was not interpreted:\n%s\n", 
            line);

  } // end input from parameter file


  /* Compute internal energies and set velocities */

  KHInnerInternalEnergy    = KHInnerPressure/((Gamma - 1.0)*KHInnerDensity);
  KHOuterInternalEnergy    = KHOuterPressure/((Gamma - 1.0)*KHOuterDensity);
  float KHInnerVelocity[3] = {0.0, 0.0, 0.0};
  float KHOuterVelocity[3] = {0.0, 0.0, 0.0};
  float KHBField[3] = {0.0, 0.0, 0.0};
  /* gas initally moving right */
  KHInnerVelocity[0]      += 0.5*KHVelocityJump + KHBulkVelocity; 
  /* gas initally moving left */
  KHOuterVelocity[0]      -= 0.5*KHVelocityJump - KHBulkVelocity; 

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  /* If KHRamp is not set, then set up ICs according to the old method.
     Two fluids separated by a discontinuity.  Initial perturbations
     are due to random fluctuations in the y-velocity of all fluid elements
     in domain.  */

  if (KHRamp == 0) {

    /* set up uniform grid without an inner flow */
    if (TopGrid.GridData->InitializeUniformGrid(KHOuterDensity, 
                                                KHOuterInternalEnergy,
                                                KHOuterInternalEnergy,
                                                KHOuterVelocity, 
                                                KHBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
    }

    /* set up the inner flow and add noise to velocities */
    if (TopGrid.GridData->KHInitializeGrid(KHInnerDensity, 
                                           KHInnerInternalEnergy,
                                           KHOuterInternalEnergy,
                                           KHPerturbationAmplitude,
                                           KHInnerVelocity[0], 
                                           KHOuterVelocity[0],
                                           KHInnerPressure,
                                           KHOuterPressure,
                                           KHRandomSeed)
        == FAIL) {
      ENZO_FAIL("Error in KHInitializeGrid.");
    }
  }

  /* If KHRamp is set, then set up ICs according to the new, ramp method.
     Two fluids separated by a continuous ramp in density and velocity.  
     Initial perturbations are due to a sinusoidal fluctuations in the 
     y-velocity of all fluid elements in the domain.  These ICs give
     convergent behavior as resolution increases.  */

  else {

    /* initialize grid and fields */
    if (TopGrid.GridData->InitializeUniformGrid(KHOuterDensity, 
                                                KHOuterInternalEnergy,
                                                KHOuterInternalEnergy,
                                                KHOuterVelocity, 
                                                KHBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
    }

    if (TopGrid.GridData->KHInitializeGridRamp(KHInnerDensity, 
                                               KHOuterDensity,
                                               KHInnerInternalEnergy,
                                               KHOuterInternalEnergy,
                                               KHPerturbationAmplitude,
                                               KHInnerVelocity[0], 
                                               KHOuterVelocity[0],
                                               KHInnerPressure,
                                               KHOuterPressure,
                                               KHRampWidth)
        == FAIL) {
      ENZO_FAIL("Error in KHInitializeGridRamp.");
    }

    /* If requested, refine the grid to the desired level. */

    if (RefineAtStart) {

      /* Declare, initialize and fill out the LevelArray. */

      LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
        LevelArray[level] = NULL;
      AddLevel(LevelArray, &TopGrid, 0);

      /* Add levels to the maximum depth or until no new levels are created,
         and re-initialize the level after it is created. */

      for (level = 0; level < MaximumRefinementLevel; level++) {
        if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
          ENZO_FAIL("Error in RebuildHierarchy.");
        }
        if (LevelArray[level+1] == NULL)
          break;
        LevelHierarchyEntry *Temp = LevelArray[level+1];
        while (Temp != NULL) {
          if (Temp->GridData->KHInitializeGridRamp(KHInnerDensity, 
                                                   KHOuterDensity,
                                                   KHInnerInternalEnergy,
                                                   KHOuterInternalEnergy,
                                                   KHPerturbationAmplitude,
                                                   KHInnerVelocity[0], 
                                                   KHOuterVelocity[0],
                                                   KHInnerPressure,
                                                   KHOuterPressure,
                                                   KHRampWidth)
              == FAIL) {
            ENZO_FAIL("Error in KHInitializeGridRamp.");
          }
          Temp = Temp->NextGridThisLevel;
        }
      } // end: loop over levels

      /* Loop back from the bottom, restoring the consistency among levels. */

      for (level = MaximumRefinementLevel; level > 0; level--) {
        LevelHierarchyEntry *Temp = LevelArray[level];
        while (Temp != NULL) {
          if (Temp->GridData->ProjectSolutionToParentGrid(
                   *LevelArray[level-1]->GridData) == FAIL) {
            ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
          }
          Temp = Temp->NextGridThisLevel;
        }
      }
    } // end: if RefineAtStart
  }
  printf("KH: single grid start-up.\n");

  /* set up field names and units */

  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "KHInnerDensity  = %"FSYM"\n", KHInnerDensity);
    fprintf(Outfptr, "KHInnerPressure = %"FSYM"\n", KHInnerPressure);
    fprintf(Outfptr, "KHOuterDensity  = %"FSYM"\n", KHOuterDensity);
    fprintf(Outfptr, "KHOuterPressure = %"FSYM"\n", KHOuterPressure);
  }

  return SUCCESS;

}
