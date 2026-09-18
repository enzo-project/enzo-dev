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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int AGNDiskInitialize(FILE *fptr, FILE *Outfptr, 
		      HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *BxName = "Bx";
  const char *ByName = "By";
  const char *BzName = "Bz";
  const char *PhiName = "Phi";
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
  const char *kphHIName    = "HI_kph";
  const char *gammaHIName  = "HI_gamma";
  const char *kphHeIName   = "HeI_kph";
  const char *gammaHeIName = "HeI_gamma";
  const char *kphHeIIName  = "HeII_kph";
  const char *gammaHeIIName= "HeII_gamma";
  const char *kdissH2IName = "H2I_kdiss";
  const char *RadAccel1Name = "RadAccel1";
  const char *RadAccel2Name = "RadAccel2";
  const char *RadAccel3Name = "RadAccel3";
  const char *GravPotenName = "PotentialField";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int level, i;

  /* set default parameters */

  int RefineAtStart = 0;
  int DiskType = 1, BlackHoleType = 0;
  float DiskDensity = 1.0,
    DiskTemperature = 1.0,
    BlackHoleMass = 0.0; 
  FLOAT DiskRadius = 1.0, 
    DiskHeight = 1.0;
  int UseGas = 1;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {


    /* read parameters */

    sscanf(line, "DiskType = %" ISYM,
		  &DiskType);
    sscanf(line, "RefineAtStart = %" ISYM, 
		  &RefineAtStart);
    sscanf(line, "BlackHoleMass = %" FSYM,
		  &BlackHoleMass);
    sscanf(line, "DiskType = %" ISYM, 
		  &DiskType);
    sscanf(line, "BlackHoleType = %" ISYM, 
		  &BlackHoleType);
    sscanf(line, "UseGas = %" ISYM, 
		  &UseGas);
    sscanf(line, "DiskDensity = %" FSYM, 
		  &DiskDensity);
    sscanf(line, "DiskTemperature = %" FSYM, 
		  &DiskTemperature);
    sscanf(line, "DiskRadius = %" PSYM, 
		  &DiskRadius);
    sscanf(line, "DiskHeight = %" PSYM, 
		  &DiskHeight);

  } // end input from parameter file

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, 1);
  float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;
  float PressureUnits = DensityUnits*pow(VelocityUnits,2);

  printf("denu=%" GSYM", velu=%" GSYM", timeu=%" GSYM", tempu=%" GSYM", lenu=%" GSYM", bu=%" GSYM", presu=%" GSYM"\n",
	 DensityUnits, VelocityUnits, TimeUnits, TemperatureUnits, LengthUnits, 
	 MagneticUnits, PressureUnits);

  printf("timeu=%" GSYM"(year)\n", TimeUnits/3.1558e7);
  //printf("temp=%" GSYM", radius=%" GSYM", height=%" GSYM", density=%" GSYM"\n",
  // DiskTemperature[0], DiskRadius[0], DiskHeight[0], DiskDensity[0]);

  if (UsePhysicalUnit) {
    DiskDensity /= DensityUnits;
    DiskTemperature /= TemperatureUnits;
  }

  /* set up grid */

  if (TopGrid.GridData->AGNDiskInitializeGrid(
	     BlackHoleMass, BlackHoleType, DiskType,
	     DiskDensity,
	     DiskTemperature,
	     DiskRadius, DiskHeight, UseGas, 0) == FAIL) {
    fprintf(stderr, "Error in GalaxyDiskInitializeGrid.\n");
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
	if (Temp->GridData->AGNDiskInitializeGrid(
	     BlackHoleMass, BlackHoleType, DiskType,
	     DiskDensity,
	     DiskTemperature,
	     DiskRadius, DiskHeight, UseGas, level+1) == FAIL) {
	  fprintf(stderr, "Error in GalaxyDiskInitializeGrid.\n");
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
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  if (HydroMethod == MHD_RK) {
      DataLabel[count++] = BxName;
      DataLabel[count++] = ByName;
      DataLabel[count++] = BzName;
      DataLabel[count++] = PhiName;
  }
  if (MultiSpecies) {
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = HMName;
      DataLabel[count++] = H2IName;
      DataLabel[count++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = DIName;
      DataLabel[count++] = DIIName;
      DataLabel[count++] = HDIName;
    }
  }  // if Multispecies                                                                                                      
#ifdef TRANSFER
  if (RadiativeTransfer) {
    if (MultiSpecies) {
      DataLabel[count++]  = kphHIName;
      DataLabel[count++]  = gammaHIName;
      DataLabel[count++]  = kphHeIName;
      DataLabel[count++]  = gammaHeIName;
      DataLabel[count++]  = kphHeIIName;
      DataLabel[count++]  = gammaHeIIName;
      if (MultiSpecies > 1)
        DataLabel[count++]= kdissH2IName;
    }
  }

  if (RadiationPressure) {
    DataLabel[count++]  = RadAccel1Name;
    DataLabel[count++]  = RadAccel2Name;
    DataLabel[count++]  = RadAccel3Name;
  }
#endif
  if (WritePotential) {
    DataLabel[count++] = GravPotenName;
  }


  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  return SUCCESS;

}
