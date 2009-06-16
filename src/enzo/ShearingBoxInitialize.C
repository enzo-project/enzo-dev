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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, EFLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, EFLOAT Time);

int ShearingBoxInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const  char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *BxName = "Bx";
  const char *ByName = "By";
  const char *BzName = "Bz";
  const char *PhiName = "Phi";
  const char *DebugName = "Debug";
  const char *Phi_pName = "Phip";
  
  /* declarations */

  char  line[MAX_LINE_LENGTH];
  

  /* set default parameters */

 
  /* read input from file */
  int ShearingBoxProblemType = 0; // 0 = advecting sphere; 1 = shearing box

  float ThermalMagneticRatio=400; 
  float FluctuationAmplitudeFraction=0.1;
  int ShearingBoxRefineAtStart   = FALSE;
  float ShearingGeometry=0.5;
  int InitialMagneticFieldConfiguration=0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "ShearingBoxRefineAtStart = %d", &RefineAtStart); 
    ret += sscanf(line, "ThermalMagneticRatio= %f", &ThermalMagneticRatio);
    ret += sscanf(line, "FluctuationAmplitudeFraction = %f", &FluctuationAmplitudeFraction);
    ret += sscanf(line, "ShearingBoxProblemType = %d", &ShearingBoxProblemType);  

  } 


  if (TopGrid.GridData->ShearingBoxInitializeGrid(ThermalMagneticRatio, FluctuationAmplitudeFraction, ShearingGeometry, ShearingBoxProblemType, InitialMagneticFieldConfiguration)
      == FAIL) 
    ENZO_FAIL("Error in ShearingBoxInitializeGrid.\n");
  


  if (ShearingBoxRefineAtStart) {
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
	ENZO_FAIL("");
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->ShearingBoxInitializeGrid(ThermalMagneticRatio, FluctuationAmplitudeFraction, ShearingGeometry, ShearingBoxProblemType, InitialMagneticFieldConfiguration)
	    == FAIL) 
	  ENZO_FAIL("Error in ShearingBoxInitializeGrid.\n");
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels
  } // end: if (CollapseTestRefineAtStart)




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
  if(useMHD){
  DataLabel[count++] = BxName;
  DataLabel[count++] = ByName;
  DataLabel[count++] = BzName;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return ENZO_SUCCESS;

}
