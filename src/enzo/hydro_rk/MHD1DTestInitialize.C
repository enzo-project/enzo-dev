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
  float  LeftDensity = 1.0, RightDensity = 1.0, 
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
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "LeftVelocityX = %"FSYM,
		  &LeftVelocityX);
    ret += sscanf(line, "LeftVelocityY = %"FSYM,
		  &LeftVelocityY);
    ret += sscanf(line, "LeftVelocityZ = %"FSYM,
		  &LeftVelocityZ);
    ret += sscanf(line, "LeftPressure = %"FSYM, 
		  &LeftPressure);
    ret += sscanf(line, "LeftDensity = %"FSYM, 
		  &LeftDensity);
    ret += sscanf(line, "LeftBx = %"FSYM,
		  &LeftBx);
    ret += sscanf(line, "LeftBy = %"FSYM,
		  &LeftBy);
    ret += sscanf(line, "LeftBz = %"FSYM,
		  &LeftBz);
    ret += sscanf(line, "RightVelocityX = %"FSYM, 
		  &RightVelocityX);
    ret += sscanf(line, "RightVelocityY = %"FSYM, 
		  &RightVelocityY);
    ret += sscanf(line, "RightVelocityZ = %"FSYM,
		  &RightVelocityZ);
    ret += sscanf(line, "RightPressure = %"FSYM, 
		  &RightPressure);
    ret += sscanf(line, "RightDensity = %"FSYM,
                  &RightDensity);
    ret += sscanf(line, "RightBx = %"FSYM,
		  &RightBx);
    ret += sscanf(line, "RightBy = %"FSYM,
		  &RightBy);
    ret += sscanf(line, "RightBz = %"FSYM,
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

  printf("DensityUnits=%"GSYM",VelocityUnits=%"GSYM",LengthUnits=%"GSYM",TimeUnits=%"GSYM" (%"GSYM" yr),PressureUnits=%"GSYM"\n", 
	 DensityUnits, VelocityUnits, LengthUnits, TimeUnits, TimeUnits/3.1558e7, PressureUnits);

  if (UsePhysicalUnit) {
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

  if (TopGrid.GridData->MHD1DTestInitializeGrid(LeftDensity, RightDensity,
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
      printf("In level %"ISYM"\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->MHD1DTestInitializeGrid(LeftDensity, RightDensity,
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
  

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RefineAtStart      = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "LeftDensity       = %"FSYM"\n",
	    LeftDensity);
    fprintf(Outfptr, "RightDensity          = %"FSYM"\n",
	    RightDensity);
    fprintf(Outfptr, "LeftVelocityX = %"FSYM"\n",
	    LeftVelocityX);
    fprintf(Outfptr, "RightVelocityX = %"FSYM"\n",
            RightVelocityX);
    fprintf(Outfptr, "LeftVelocityY = %"FSYM"\n",
	    LeftVelocityY);
    fprintf(Outfptr, "RightVelocityY = %"FSYM"\n",
            RightVelocityY);
    fprintf(Outfptr, "LeftBx = %"FSYM"\n",
	    LeftBx);
    fprintf(Outfptr, "RightBx = %"FSYM"\n",
	    RightBx);
    fprintf(Outfptr, "LeftBy = %"FSYM"\n",
	    LeftBy);
    fprintf(Outfptr, "RightBy = %"FSYM"\n",
	    RightBy);
    fprintf(Outfptr, "LeftPressure = %"FSYM"\n",
            LeftPressure);
    fprintf(Outfptr, "RightPressure = %"FSYM"\n",
            RightPressure);
  }
  //return FAIL;
  return SUCCESS;

}

