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
  char *DensName = "Density";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *kphHIName    = "HI_kph";
  char *gammaHIName  = "HI_gamma";
  char *kphHeIName   = "HeI_kph";
  char *gammaHeIName = "HeI_gamma";
  char *kphHeIIName  = "HeII_kph";
  char *gammaHeIIName= "HeII_gamma";
  char *kdissH2IName = "H2I_kdiss";
  char *RadAccel1Name = "RadAccel1";
  char *RadAccel2Name = "RadAccel2";
  char *RadAccel3Name = "RadAccel3";
  char *GravPotenName = "PotentialField";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

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

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "DiskType = %"ISYM,
		  &DiskType);
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "BlackHoleMass = %"FSYM,
		  &BlackHoleMass);
    ret += sscanf(line, "DiskType = %"ISYM, 
		  &DiskType);
    ret += sscanf(line, "BlackHoleType = %"ISYM, 
		  &BlackHoleType);
    ret += sscanf(line, "UseGas = %"ISYM, 
		  &UseGas);
    ret += sscanf(line, "DiskDensity = %"FSYM, 
		  &DiskDensity);
    ret += sscanf(line, "DiskTemperature = %"FSYM, 
		  &DiskTemperature);
    ret += sscanf(line, "DiskRadius = %"PSYM, 
		  &DiskRadius);
    ret += sscanf(line, "DiskHeight = %"PSYM, 
		  &DiskHeight);

  } // end input from parameter file

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, 1);
  float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;
  float PressureUnits = DensityUnits*pow(VelocityUnits,2);

  printf("denu=%"GSYM", velu=%"GSYM", timeu=%"GSYM", tempu=%"GSYM", lenu=%"GSYM", bu=%"GSYM", presu=%"GSYM"\n",
	 DensityUnits, VelocityUnits, TimeUnits, TemperatureUnits, LengthUnits, 
	 MagneticUnits, PressureUnits);

  printf("timeu=%"GSYM"(year)\n", TimeUnits/3.1558e7);
  //printf("temp=%"GSYM", radius=%"GSYM", height=%"GSYM", density=%"GSYM"\n",
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
