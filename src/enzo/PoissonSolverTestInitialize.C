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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int PoissonSolverTestInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  if (!UseMHD) 
    ENZO_FAIL("DivergenceCleaning only useful with MHD simulations");

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
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int TestType=0;
  float TestGeometryControl=1;
  int RefineAtStart=0;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "PoissonSolverTestType = %d", &TestType);
    ret += sscanf(line, "PoissonSolverTestGeometryControl = %f", &TestGeometryControl);
    ret += sscanf(line, "PoissonSolverTestRefineAtStart = %d", &RefineAtStart);
    /* if the line is suspicious, issue a warning */

  } // end input from parameter file
  
 

  if (TopGrid.GridData->PoissonSolverTestInitializeGrid(TestType, TestGeometryControl) == FAIL) {
    fprintf(stderr, "Error in PoissonSolverTestInitializeGrid.\n");
    return FAIL;
  }

  /* If requested, refine the grid to the desired level. */


  if (RefineAtStart) {
     /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (int level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	ENZO_FAIL("Error in RebuildHierarchy.\n");
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->PoissonSolverTestInitializeGrid(TestType, TestGeometryControl)
	    == FAIL) 
	  ENZO_FAIL("Error in ShearingBoxInitializeGrid.\n");
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels
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
  DataLabel[count++] = BxName;
  DataLabel[count++] = ByName;
  DataLabel[count++] = BzName;
  DataLabel[count++] = PhiName;
  if(UseDivergenceCleaning){

    DataLabel[count++] = Phi_pName;
  }


  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return SUCCESS;

}
