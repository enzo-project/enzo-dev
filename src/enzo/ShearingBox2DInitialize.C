/***********************************************************************
/
/  INITIALIZE A SHEARING BOX TEST
/
/  written by: Fen Zhao
/  date:       June, 2009
/  modified1:
/
/  PURPOSE:
/    Set up exither an advecting sphere or the standard shearing box simluation
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/



#include <stdlib.h>
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

// void WriteListOfFloats(FILE *fptr, int N, float floats[]);
// void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
// void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
// int RebuildHierarchy(TopGridData *MetaData,
// 		     LevelHierarchyEntry *LevelArray[], int level);
// void WriteListOfFloats(FILE *fptr, int N, float floats[]);
// void WriteListOfFloats(FILE *fptr, int N, EFLOAT floats[]);
// void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
// int RebuildHierarchy(TopGridData *MetaData,
// 		     LevelHierarchyEntry *LevelArray[], int level);
// int GetUnits(float *DensityUnits, float *LengthUnits,
// 		      float *TemperatureUnits, float *TimeUnits,
// 		      float *VelocityUnits, EFLOAT Time);

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);


int ShearingBox2DInitialize (FILE *fptr, FILE *Outfptr, 
			       HierarchyEntry &TopGrid, TopGridData &MetaData)

{

  //Lets Assume the shear is in the y direction and the x direction is radial

  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const  char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  //const char *Vel3Name = "z-velocity";
  const char *BxName = "Bx";
  const char *ByName = "By";
  // const char *BzName = "Bz";
  const char *PhiName = "Phi";
  const char *DebugName = "Debug";
  const char *Phi_pName = "Phip";
  
  /* declarations */

  char  line[MAX_LINE_LENGTH];
  

  /* set default parameters */

 
  /* read input from file */
 

  float ThermalMagneticRatio=400; 
  float FluctuationAmplitudeFraction=0.1;
  int ShearingBoxRefineAtStart   = FALSE;
  float ShearingGeometry=2.0;
  int InitialMagneticFieldConfiguration=0;
  int RefineAtStart=1;

  int ret;

  

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

/* read parameters */
    ret += sscanf(line, "ShearingBoxRefineAtStart = %"ISYM, &RefineAtStart); 
    ret += sscanf(line, "ShearingBoxThermalMagneticRatio= %"FSYM, &ThermalMagneticRatio);
    ret += sscanf(line, "ShearingBoxFluctuationAmplitudeFraction = %"FSYM, &FluctuationAmplitudeFraction);
    ret += sscanf(line, "ShearingBoxGeometry = %"FSYM, &ShearingGeometry);  
    ret += sscanf(line, "ShearingBoxInitialMagneticFieldConfiguration = %"ISYM, &InitialMagneticFieldConfiguration);  
 
  } 
 


  if (TopGrid.GridData->ShearingBox2DInitializeGrid(ThermalMagneticRatio, FluctuationAmplitudeFraction, ShearingGeometry, InitialMagneticFieldConfiguration)
      == FAIL) 
    ENZO_FAIL("Error in ShearingBoxInitializeGrid.\n");
  


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
	if (Temp->GridData->ShearingBox2DInitializeGrid(ThermalMagneticRatio, FluctuationAmplitudeFraction, ShearingGeometry, InitialMagneticFieldConfiguration)
	    == FAIL) 
	  ENZO_FAIL("Error in ShearingBoxInitializeGrid.\n");
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels
  } // end: if (CollapseTestRefineAtStart)




  /* set up field names and units */

  int count = 0;
  DataLabel[count++] =  (char*) DensName;
  DataLabel[count++] =  (char*) Vel1Name;
  DataLabel[count++] =  (char*) Vel2Name;
  //DataLabel[count++] =  (char*) Vel3Name;
  DataLabel[count++] =  (char*) TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] =  (char*) GEName;
  }
  if(UseMHD){

  DataLabel[count++] =  (char*) BxName;
  DataLabel[count++] =  (char*) ByName;
  // DataLabel[count++] =  (char*) BzName;
  DataLabel[count++] =  (char*) PhiName;
  }

  for (int i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }




  return SUCCESS;

}
