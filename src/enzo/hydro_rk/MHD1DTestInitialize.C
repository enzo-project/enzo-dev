/***********************************************************************
/
/  INITIALIZE MHD 1D TEST
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
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

int MHD1DTestInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid, TopGridData &MetaData) 
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

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  float  RampWidth,
    LeftDensity = 1.0, RightDensity = 1.0, 
    LeftVelocityX = 0, RightVelocityX = 0, 
    LeftVelocityY = 0, RightVelocityY = 0, 
    LeftVelocityZ = 0, RightVelocityZ = 0,
    LeftPressure = 1.0, RightPressure = 1.0,
    LeftBx = 0.0, RightBx = 0.0,
    LeftBy = 0.0, RightBy = 0.0,
    LeftBz = 0.0, RightBz = 0.0;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "RampWidth = %f",
		  &RampWidth);
    ret += sscanf(line, "RefineAtStart = %d", 
		  &RefineAtStart);
    ret += sscanf(line, "LeftVelocityX = %f",
		  &LeftVelocityX);
    ret += sscanf(line, "LeftVelocityY = %f",
		  &LeftVelocityY);
    ret += sscanf(line, "LeftVelocityZ = %f",
		  &LeftVelocityZ);
    ret += sscanf(line, "LeftPressure = %f", 
		  &LeftPressure);
    ret += sscanf(line, "LeftDensity = %f", 
		  &LeftDensity);
    ret += sscanf(line, "LeftBx = %f",
		  &LeftBx);
    ret += sscanf(line, "LeftBy = %f",
		  &LeftBy);
    ret += sscanf(line, "LeftBz = %f",
		  &LeftBz);
    ret += sscanf(line, "RightVelocityX = %f", 
		  &RightVelocityX);
    ret += sscanf(line, "RightVelocityY = %f", 
		  &RightVelocityY);
    ret += sscanf(line, "RightVelocityZ = %f",
		  &RightVelocityZ);
    ret += sscanf(line, "RightPressure = %f", 
		  &RightPressure);
    ret += sscanf(line, "RightDensity = %f",
                  &RightDensity);
    ret += sscanf(line, "RightBx = %f",
		  &RightBx);
    ret += sscanf(line, "RightBy = %f",
		  &RightBy);
    ret += sscanf(line, "RightBz = %f",
		  &RightBz);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
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
  double PressureUnits = DensityUnits*pow(VelocityUnits,2);
  double MagneticUnits = sqrt(4.0*M_PI*PressureUnits);

  printf("DensityUnits=%g,VelocityUnits=%g,LengthUnits=%g,TimeUnits=%g (%g yr),PressureUnits=%g\n", 
	 DensityUnits, VelocityUnits, LengthUnits, TimeUnits, TimeUnits/3.1558e7, PressureUnits);

  if (UsePhysicalUnit) {
    RampWidth /= LengthUnits;
    LeftDensity /= DensityUnits;
    RightDensity /= DensityUnits;
    LeftVelocityX  /= VelocityUnits;
    RightVelocityX  /= VelocityUnits;
    LeftVelocityY  /= VelocityUnits;
    RightVelocityY  /= VelocityUnits;
    LeftVelocityZ  /= VelocityUnits;
    RightVelocityZ  /= VelocityUnits;
    LeftPressure   /= PressureUnits;
    RightPressure   /= PressureUnits;
    LeftBx  /= MagneticUnits;
    RightBx  /= MagneticUnits;
    LeftBy  /= MagneticUnits;
    RightBy  /= MagneticUnits;
    LeftBz  /= MagneticUnits;
    RightBz  /= MagneticUnits;
  }

  /* set up grid */

  if (TopGrid.GridData->MHD1DTestInitializeGrid(RampWidth,
						LeftDensity, RightDensity,
						LeftVelocityX,  RightVelocityX,
						LeftVelocityY,  RightVelocityY,
						LeftVelocityZ,  RightVelocityZ,
						LeftPressure,   RightPressure,
						LeftBx,  RightBx,
						LeftBy,  RightBy,
						LeftBz,  RightBz)  == FAIL) {
    fprintf(stderr, "Error in MHD1DTestInitializeGrid.\n");
    return FAIL;
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

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
      printf("In level %i\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->MHD1DTestInitializeGrid(RampWidth,
						    LeftDensity, RightDensity,
						    LeftVelocityX,  RightVelocityX,
						    LeftVelocityY,  RightVelocityY,
						    LeftVelocityZ,  RightVelocityZ, 
						    LeftPressure,   RightPressure,
						    LeftBx,  RightBx,
						    LeftBy,  RightBy,
						    LeftBz,  RightBz) == FAIL) {
	  fprintf(stderr, "Error in MHD1DTestInitializeGrid.\n");
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
  DataLabel[count++] = BxName;
  DataLabel[count++] = ByName;
  DataLabel[count++] = BzName;
  DataLabel[count++] = PhiName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RefineAtStart      = %d\n",
	    RefineAtStart);
    fprintf(Outfptr, "RampWidth         = %f\n",
	    RampWidth);
    fprintf(Outfptr, "LeftDensity       = %f\n",
	    LeftDensity);
    fprintf(Outfptr, "RightDensity          = %f\n",
	    RightDensity);
    fprintf(Outfptr, "LeftVelocityX = %f\n",
	    LeftVelocityX);
    fprintf(Outfptr, "RightVelocityX = %f\n",
            RightVelocityX);
    fprintf(Outfptr, "LeftVelocityY = %f\n",
	    LeftVelocityY);
    fprintf(Outfptr, "RightVelocityY = %f\n",
            RightVelocityY);
    fprintf(Outfptr, "LeftBx = %f\n",
	    LeftBx);
    fprintf(Outfptr, "RightBx = %f\n",
	    RightBx);
    fprintf(Outfptr, "LeftBy = %f\n",
	    LeftBy);
    fprintf(Outfptr, "RightBy = %f\n",
	    RightBy);
    fprintf(Outfptr, "LeftPressure = %f\n",
            LeftPressure);
    fprintf(Outfptr, "RightPressure = %f\n",
            RightPressure);
  }
  //return FAIL;
  return SUCCESS;

}

