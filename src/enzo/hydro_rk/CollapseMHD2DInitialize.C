/***********************************************************************
/
/  INITIALIZE MAGNETIZED CLOUD
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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CollapseMHD2DInitialize(FILE *fptr, FILE *Outfptr, 
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
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int n_sphere = 1;
  int RefineAtStart   = TRUE;
  int UseParticles    = FALSE;
  float MediumDensity = 1.0, 
    MediumPressure = 1.0;
  int   SphereType;
  float SphereDensity = 1.0,
    SpherePressure = 1.0,
    SphereSoundVelocity = 1.0,
    SphereAngVel = 0.0,
    SphereCutOff = 1.0;
  FLOAT SphereRadius = 1.0, SphereCoreRadius = 1.0;
  float B0 = 0.0;
  float theta_B = 0.0;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "RefineAtStart = %"ISYM, &RefineAtStart);
    ret += sscanf(line, "MediumDensity = %"FSYM, &MediumDensity);
    ret += sscanf(line, "MediumPressure = %"FSYM, &MediumPressure);
    ret += sscanf(line, "InitialBField = %"FSYM, &B0);
    ret += sscanf(line, "theta_B = %"FSYM, &theta_B);
 
    ret += sscanf(line, "SphereType = %"ISYM, &SphereType);
    ret += sscanf(line, "SphereRadius = %"PSYM, &SphereRadius);
    ret += sscanf(line, "SphereCoreRadius = %"PSYM, &SphereCoreRadius);
    ret += sscanf(line, "SphereDensity = %"FSYM, &SphereDensity);
    ret += sscanf(line, "SpherePressure = %"FSYM, &SpherePressure);
    ret += sscanf(line, "SphereSoundVelocity = %"FSYM, &SphereSoundVelocity);
    ret += sscanf(line, "SphereAngVel = %"FSYM, &SphereAngVel);
    ret += sscanf(line, "SphereCutOff = %"FSYM, &SphereCutOff);

  } // end input from parameter file
  
  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0, bfieldu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;
    bfieldu = sqrt(presu*4.0*Pi);
  }
  
  printf("rhou=%"GSYM",velu=%"GSYM",lenu=%"GSYM",tu=%"GSYM" (%"GSYM" yr),tempu=%"GSYM",presu=%"GSYM", bfieldu=%"GSYM", tempu=%"GSYM"\n", 
	 rhou, velu,lenu,tu,tu/3.1558e7,tempu,presu,bfieldu, tempu);

  MediumDensity /= rhou;
  MediumPressure /= presu;

  B0 /= bfieldu;

  SphereDensity /= rhou;
  SpherePressure /= presu;
  SphereSoundVelocity /= velu;
  SphereAngVel *= tu;

  printf("rhoc=%"GSYM", rhom=%"GSYM", pm=%"GSYM"\n", SphereDensity, MediumDensity, MediumPressure);


  if (TopGrid.GridData->CollapseMHD3DInitializeGrid(
	     SphereRadius, SphereCoreRadius, SphereDensity,
	     SpherePressure, SphereSoundVelocity, SphereAngVel, B0, theta_B,
	     SphereType, MediumDensity, MediumPressure, 0) == FAIL) {
    fprintf(stderr, "Error in CollapseTestInitializeGrid.\n");
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

      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->CollapseMHD3DInitializeGrid(
				SphereRadius, SphereCoreRadius, SphereDensity,
				SpherePressure, SphereSoundVelocity, SphereAngVel, B0,theta_B,
				SphereType, MediumDensity, MediumPressure, level+1) == FAIL) {
	  fprintf(stderr, "Error in Collapse3DInitializeGrid.\n");
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
    DataLabel[count++] = DebugName;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return SUCCESS;

}
