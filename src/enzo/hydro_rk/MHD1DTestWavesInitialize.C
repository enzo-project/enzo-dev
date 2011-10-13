/***********************************************************************
/
/  INITIALIZE MHD WAVE 1D TEST
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1: J. S. Oishi
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

int MHD1DTestWavesInitialize(FILE *fptr, FILE *Outfptr, 
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
  float  DensityPert = 1.0e-6, 
    VelocityXPert = 0, 
    VelocityYPert = 0, 
    VelocityZPert = 0,
    EtotPert = 1.0,
    BxPert = 0.0,
    ByPert = 0.0,
    BzPert = 0.0;
  
  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "VelocityXPert = %"FSYM,
		  &VelocityXPert);
    ret += sscanf(line, "VelocityYPert = %"FSYM,
		  &VelocityYPert);
    ret += sscanf(line, "VelocityZPert = %"FSYM,
		  &VelocityZPert);
    ret += sscanf(line, "EtotPert = %"FSYM, 
		  &EtotPert);
    ret += sscanf(line, "DensityPert = %"FSYM, 
		  &DensityPert);
    ret += sscanf(line, "BxPert = %"FSYM,
		  &BxPert);
    ret += sscanf(line, "ByPert = %"FSYM,
		  &ByPert);
    ret += sscanf(line, "BzPert = %"FSYM,
		  &BzPert);
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
    DensityPert /= DensityUnits;
    VelocityXPert  /= VelocityUnits;
    VelocityYPert  /= VelocityUnits;
    VelocityZPert  /= VelocityUnits;
    EtotPert   /= PressureUnits;
    BxPert  /= MagneticUnits;
    ByPert  /= MagneticUnits;
    BzPert  /= MagneticUnits;
  }

  /* set up grid */

  if (TopGrid.GridData->MHD1DTestWavesInitializeGrid(DensityPert, 
						VelocityXPert,
						VelocityYPert,
						VelocityZPert,
						EtotPert, 
						BxPert,
						ByPert,
						BzPert)  == FAIL) {
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
	if (Temp->GridData->MHD1DTestWavesInitializeGrid(DensityPert,    
						    VelocityXPert,  
						    VelocityYPert,  
						    VelocityZPert,  
						    EtotPert,
						    BxPert,
						    ByPert,
						    BzPert) == FAIL) {
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
    fprintf(Outfptr, "RefineAtStart      = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "DensityPert       = %"FSYM"\n",
	    DensityPert);
    fprintf(Outfptr, "VelocityXPert = %"FSYM"\n",
	    VelocityXPert);
    fprintf(Outfptr, "VelocityYPert = %"FSYM"\n",
	    VelocityYPert);
    fprintf(Outfptr, "VelocityZPert = %"FSYM"\n",
	    VelocityZPert);
    fprintf(Outfptr, "BxPert = %"FSYM"\n",
	    BxPert);
    fprintf(Outfptr, "ByPert = %"FSYM"\n",
	    ByPert);
    fprintf(Outfptr, "BzPert = %"FSYM"\n",
	    BzPert);
    fprintf(Outfptr, "EtotPert = %"FSYM"\n",
            EtotPert);
  }
  //return FAIL;
  return SUCCESS;

}

