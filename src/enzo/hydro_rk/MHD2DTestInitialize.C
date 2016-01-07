/***********************************************************************
/
/  INITIALIZE MHD 2D TEST
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1: Tom Abel 2010 
/            added many new tests including the Wengen Coliding Flow test
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
#include "LevelHierarchy.h"
#include "TopGridData.h"
#define DEFINE_STORAGE
#include "MHD2DTestGlobalData.h"
#undef DEFINE_STORAGE

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

int MHD2DTestInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid,
			TopGridData &MetaData, int SetBaryonFields) 
{
  char *DensName = "Density";
  char *PresName = "Pressure";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  int MHD2DProblemType = 0;
  RampWidth = 0.05;
  LowerDensity = 1.0; UpperDensity = 1.0;
  LowerVelocityX = 0; UpperVelocityX = 0;
  LowerVelocityY = 0; UpperVelocityY = 0; 
  LowerPressure = 1.0; UpperPressure = 1.0;
  LowerBx = 0.0; UpperBx = 0.0;
  LowerBy = 0.0; UpperBy = 0.0;
  UseColour = FALSE;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "LowerVelocityX = %"FSYM,
		  &LowerVelocityX);
    ret += sscanf(line, "LowerVelocityY = %"FSYM,
		  &LowerVelocityY);
    ret += sscanf(line, "LowerPressure = %"FSYM, 
		  &LowerPressure);
    ret += sscanf(line, "LowerDensity = %"FSYM, 
		  &LowerDensity);
    ret += sscanf(line, "LowerBx = %"FSYM,
		  &LowerBx);
    ret += sscanf(line, "LowerBy = %"FSYM,
		  &LowerBy);
    ret += sscanf(line, "UpperVelocityX = %"FSYM, 
		  &UpperVelocityX);
    ret += sscanf(line, "UpperVelocityY = %"FSYM, 
		  &UpperVelocityY);
    ret += sscanf(line, "UpperPressure = %"FSYM, 
		  &UpperPressure);
    ret += sscanf(line, "UpperDensity = %"FSYM,
                  &UpperDensity);
    ret += sscanf(line, "UpperBx = %"FSYM,
		  &UpperBx);
    ret += sscanf(line, "UpperBy = %"FSYM,
		  &UpperBy);
    ret += sscanf(line, "MHD2DProblemType = %"ISYM,
		  &MHD2DProblemType);
    ret += sscanf(line, "RampWidth = %"FSYM,
		  &RampWidth);
    ret += sscanf(line, "UseColour = %"ISYM, 
		  &UseColour);
    
    //        fprintf(stderr, "%"ISYM" MHD2DTestInitialize !!!!!!!!!!\n", RefineAtStart);
    /* if the line is suspicious, issue a warning */

    /*
    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
    */

  } // end input from parameter file



  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set up grid */
  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  int count = 0;
  while (CurrentGrid != NULL) {
    printf("count %i %i\n", count++, MyProcessorNumber);
    if (CurrentGrid->GridData->MHD2DTestInitializeGrid(MHD2DProblemType, UseColour,
						RampWidth,
						LowerDensity, UpperDensity,
						LowerVelocityX,  UpperVelocityX,
						LowerVelocityY,  UpperVelocityY,
						LowerPressure,   UpperPressure,
						LowerBx,  UpperBx,
						LowerBy,  UpperBy,
						SetBaryonFields)  == FAIL) {
    fprintf(stderr, "Error in MHD2DTestInitializeGrid.\n");
    return FAIL;
    }
    CurrentGrid = CurrentGrid->NextGridThisLevel;
    
  } // while CurrentGrid

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */


  if (SetBaryonFields && RefineAtStart) {

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
	if (Temp->GridData->MHD2DTestInitializeGrid(MHD2DProblemType, UseColour,
						    RampWidth,
						    LowerDensity, UpperDensity,
						    LowerVelocityX,  UpperVelocityX,
						    LowerVelocityY,  UpperVelocityY,
						    LowerPressure,   UpperPressure,
						    LowerBx,  UpperBx,
						    LowerBy,  UpperBy,
						    SetBaryonFields) == FAIL) {
	  fprintf(stderr, "Error in MHD2DTestInitializeGrid.\n");
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

    //WriteAllData(MetaData.DataDumpName, MetaDaGrita.DataDumpNumber,
    //       &TopGrid, MetaData, Exterior, -1);

  } // end: if (RefineAtStart)


  /* set up field names and units */

  count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  if (UseMHD) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
    if( HydroMethod == MHD_RK ){
        DataLabel[count++] = PhiName;
    }
    if(UseDivergenceCleaning){
      DataLabel[count++] = Phi_pName;
      DataLabel[count++] = DebugName;
    }
  }
  if ( UseMHDCT ){
    MHDLabel[0] = "BxF";
    MHDLabel[1] = "ByF";
    MHDLabel[2] = "BzF";
    
    MHDeLabel[0] = "Ex";
    MHDeLabel[1] = "Ey";
    MHDeLabel[2] = "Ez";
    
    MHDUnits[0] = "None";
    MHDUnits[1] = "None";
    MHDUnits[2] = "None";
    
    MHDeUnits[0] = "None";
    MHDeUnits[1] = "None";
    MHDeUnits[2] = "None";
  }
  if (UseColour == TRUE)
    DataLabel[count++] = ColourName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  return SUCCESS;

}

