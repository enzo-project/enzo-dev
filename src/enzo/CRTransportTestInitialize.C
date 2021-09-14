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

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData, 
		 ExternalBoundary *Exterior, FLOAT WriteTime);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int CRTransportTestInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData) 
{
  char *DensName = "Density";
  char *PresName = "Pressure";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *CRName = "CREnergyDensity";
  char *ColourName = "colour";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  int TestType  = 0; 

  float Center = 0.5, GasDensity = 1.0, VelocityX = 0.0, VelocityY = 0.0, 
    VelocityZ = 0.0, GasPressure = 0.0, CREnergyDensity = 1.0,
    Bx = 0.0, By = 0.0, Bz = 0.0; 

  /* read input from file */
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;
    /* read parameters */
    
    ret += sscanf(line, "CRTransportTestType = %"ISYM, 
		  &TestType);
    ret += sscanf(line, "CRTransportTestRefineAtStart = %"ISYM, 
		  &RefineAtStart);    
    ret += sscanf(line, "CRTransportTestCenter = %"FSYM, 
		  &Center);
    ret += sscanf(line, "CRTransportTestVelocityX = %"FSYM,
		  &VelocityX);
    ret += sscanf(line, "CRTransportTestVelocityY = %"FSYM,
                  &VelocityY);
    ret += sscanf(line, "CRTransportTestVelocityZ = %"FSYM,
		                &VelocityZ);
    ret += sscanf(line, "CRTransportTestGasPressure = %"FSYM, 
		  &GasPressure);
    ret += sscanf(line, "CRTransportTestDensity = %"FSYM, 
	      	  &GasDensity);
    ret += sscanf(line, "CRTransportTestCREnergyDensity = %"FSYM,
		  &CREnergyDensity);
    ret += sscanf(line, "CRTransportTestBx = %"FSYM,
                  &Bx);
    ret += sscanf(line, "CRTransportTestBy = %"FSYM,
                  &By);
    ret += sscanf(line, "CRTransportTestBz = %"FSYM,
                  &Bz);
    
    if (ret == 0 && strstr(line, "=") && strstr(line, "CRTransportTest") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file
  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set up grid */
  TopGrid.GridData->
    CRTransportTestInitializeGrid(TestType, Center, GasDensity, VelocityX, VelocityY, VelocityZ,
			       GasPressure, CREnergyDensity, Bx, By, Bz);

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */
  printf("CRTransportTestInitialize\n");
  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
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
      printf("In level %"ISYM"\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	Temp->GridData->
	  CRTransportTestInitializeGrid(TestType, Center, GasDensity, VelocityX, 
	       VelocityY, VelocityZ, GasPressure, CREnergyDensity, Bx, By, Bz);
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

    //WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
    //       &TopGrid, MetaData, Exterior, -1);

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
    DataLabel[count++] =  BxName;
    DataLabel[count++] =  ByName;
    DataLabel[count++] =  BzName;
    DataLabel[count++] =  PhiName;
  }
  DataLabel[count++] = CRName;
  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CRTransportTestRefineAtStart        = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "CRTransportTestCenter = %"FSYM"\n",
	    Center);
    fprintf(Outfptr, "CRTransportTestDensity          = %"FSYM"\n",
	    GasDensity);
    fprintf(Outfptr, "CRTransportTestVelocityX        = %"FSYM"\n",
	    VelocityX);
    fprintf(Outfptr, "CRTransportTestVelocityY        = %"FSYM"\n",
	    VelocityY);
    fprintf(Outfptr, "CRTransportTestVelocityZ       = %"FSYM"\n",
	    VelocityZ);
    fprintf(Outfptr, "CRTransportTestGasPressure         = %"FSYM"\n",
	    GasPressure);
    fprintf(Outfptr, "CRTransportCREnergyDensity       = %"FSYM"\n",
	    CREnergyDensity);
    fprintf(Outfptr, "CRTransportTestBx     = %"FSYM"\n",
	    Bx);
    fprintf(Outfptr, "CRTransportTestBy     = %"FSYM"\n",
            By);
    fprintf(Outfptr, "CRTransportTestBz     = %"FSYM"\n",
            Bz);
  }

  return SUCCESS;

}

