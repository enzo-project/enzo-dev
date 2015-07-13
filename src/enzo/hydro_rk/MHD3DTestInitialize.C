/***********************************************************************
/
/  INITIALIZE MHD 3D TEST
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

int MHD3DTestInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid, TopGridData &MetaData) 
{
  const char *DensName = "Density";
  const char *PresName = "Pressure";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "colour";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *BxName = "Bx";
  const char *ByName = "By";
  const char *BzName = "Bz";
  const char *PhiName = "Phi";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = FALSE;
  int MHD3DProblemType = 0;
  float  rhol = 1.0, rhou = 1.0, 
    vxl = 0, vxu = 0, 
    vyl = 0, vyu = 0, 
    pl = 1.0, pu = 1.0,
    Bxl = 0.0, Bxu = 0.0,
    Byl = 0.0, Byu = 0.0;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "LeftVelocityX = %"FSYM,
		  &vxl);
    ret += sscanf(line, "LeftVelocityY = %"FSYM,
		  &vyl);
    ret += sscanf(line, "LeftPressure = %"FSYM, 
		  &pl);
    ret += sscanf(line, "LeftDensity = %"FSYM, 
		  &rhol);
    ret += sscanf(line, "LeftBx = %"FSYM,
		  &Bxl);
    ret += sscanf(line, "LeftBy = %"FSYM,
		  &Byl);
    ret += sscanf(line, "RightVelocityX = %"FSYM, 
		  &vxu);
    ret += sscanf(line, "RightVelocityY = %"FSYM, 
		  &vyu);
    ret += sscanf(line, "RightPressure = %"FSYM, 
		  &pu);
    ret += sscanf(line, "RightDensity = %"FSYM,
                  &rhou);
    ret += sscanf(line, "RightBx = %"FSYM,
		  &Bxu);
    ret += sscanf(line, "RightBy = %"FSYM,
		  &Byu);
    ret += sscanf(line, "MHD3DProblemType = %"ISYM,
		  &MHD3DProblemType);

  } // end input from parameter file

  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* set up grid */

  if (TopGrid.GridData->MHD3DTestInitializeGrid(MHD3DProblemType,
						rhol, rhou,
						vxl,  vxu,
						vyl,  vyu,
						pl,   pu,
						Bxl,  Bxu,
						Byl,  Byu)  == FAIL) {
    fprintf(stderr, "Error in MHD3DTestInitializeGrid.\n");
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
	if (Temp->GridData->MHD3DTestInitializeGrid(MHD3DProblemType,
						    rhol, rhou,
						    vxl,  vxu,
						    vyl,  vyu,
						    pl,   pu,
						    Bxl,  Bxu,
						    Byl,  Byu) == FAIL) {
	  fprintf(stderr, "Error in MHD3DTestInitializeGrid.\n");
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
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = (char*) GEName;
  }

  if (HydroMethod == MHD_RK) {
    DataLabel[count++] = (char*) BxName;
    DataLabel[count++] = (char*) ByName;
    DataLabel[count++] = (char*) BzName;
    DataLabel[count++] = (char*) PhiName;
  }

  if (MultiSpecies) {
    DataLabel[count++] = (char*) ElectronName;
    DataLabel[count++] = (char*) HIName;
    DataLabel[count++] = (char*) HIIName;
    DataLabel[count++] = (char*) HeIName;
    DataLabel[count++] = (char*) HeIIName;
    DataLabel[count++] = (char*) HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = (char*) HMName;
      DataLabel[count++] = (char*) H2IName;
      DataLabel[count++] = (char*) H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = (char*) DIName;
      DataLabel[count++] = (char*) DIIName;
      DataLabel[count++] = (char*) HDIName;
    }
  }  // if Multispecies

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  /*if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "RefineAtStart      = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "LeftDensity       = %"FSYM"\n",
	    rhol);
    fprintf(Outfptr, "RightDensity          = %"FSYM"\n",
	    rhor);
    fprintf(Outfptr, "LeftVelocityX = %"FSYM"\n",
	    vxl);
    fprintf(Outfptr, "RightVelocityX = %"FSYM"\n",
            vxr);
    fprintf(Outfptr, "LeftVelocityY = %"FSYM"\n",
	    vyl);
    fprintf(Outfptr, "RightVelocityY = %"FSYM"\n",
            vyr);
    fprintf(Outfptr, "LeftBx = %"FSYM"\n",
	    Bxl);
    fprintf(Outfptr, "RightBx = %"FSYM"\n",
	    Bxr);
    fprintf(Outfptr, "LeftBy = %"FSYM"\n",
	    Byl);
    fprintf(Outfptr, "RightBy = %"FSYM"\n",
	    Byr);
    fprintf(Outfptr, "LeftPressure = %"FSYM"\n",
            pl);
    fprintf(Outfptr, "RightPressure = %"FSYM"\n",
            pr);
	    }*/

  return SUCCESS;

}

