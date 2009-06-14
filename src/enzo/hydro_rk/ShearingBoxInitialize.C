/***********************************************************************
/
/  INITIALIZE SHEARING BOX
/
/  written by: Fen Zhao
/  date:       2008
/  modified1: Peng Wang
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
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

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int ShearingBoxInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *Phi_pName = "Phip";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int ret, level;

  /* set default parameters */

  int RefineAtStart   = TRUE;
  float AngularVelocity = 1e-3;
  float VelocityGradient = 1.5;
  float ThermalMagneticRatio = 400; 
  int ShearingBoxProblemType = 0;
 
  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "RefineAtStart = %d", &RefineAtStart);
    ret += sscanf(line, "AngularVelocity = %f", &AngularVelocity);
    ret += sscanf(line, "VelocityGradient = %f", &VelocityGradient);
    ret += sscanf(line, "ThermalMagneticRatio = %f", &ThermalMagneticRatio);
    ret += sscanf(line, "ShearingBoxProblemType = %d", &ShearingBoxProblemType);   

  } // end input from parameter file
    
  if (TopGrid.GridData->ShearingBoxInitializeGrid(AngularVelocity, VelocityGradient, 
			 ThermalMagneticRatio, ShearingBoxProblemType) == FAIL) {
    fprintf(stderr, "Error in ShearingBoxInitializeGrid.\n");
    return FAIL;
  }

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
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->ShearingBoxInitializeGrid(AngularVelocity, VelocityGradient,
			      ThermalMagneticRatio, ShearingBoxProblemType) == FAIL) {
	  fprintf(stderr, "Error in ShearingBoxInitializeGrid.\n");
	  return FAIL;
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
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (RefineAtStart)


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  if (HydroMethod == MHD_RK) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
    DataLabel[count++] = PhiName;
  }
  if(UseDivergenceCleaning){
    DataLabel[count++] = Phi_pName;
  }

  for (int i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return SUCCESS;

}
