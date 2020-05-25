/***********************************************************************
/
/  INITIALIZE ISOLATED DISK GALAXY
/
/  written by: Peng Wang
/  date:       September, 2007
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int GalaxyDiskInitialize(FILE *fptr, FILE *Outfptr, 
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
  char *Phi_pName = "Phip";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int NumberOfHalos = 1;
  int RefineAtStart   = TRUE;
  int UseParticles    = FALSE,
    UseGas = TRUE;
  float MediumTemperature = 1000, MediumDensity = 1.0;
  int   GalaxyType[MAX_SPHERES];
  float HaloDensity[MAX_SPHERES],
    HaloTemperature[MAX_SPHERES],
    HaloVelocity[MAX_SPHERES][MAX_DIMENSION],
    HaloSpin[MAX_SPHERES],
    HaloAngVel[MAX_SPHERES],
    DiskDensity[MAX_SPHERES],
    DiskTemperature[MAX_SPHERES],
    DiskMassFraction[MAX_SPHERES],
    DiskFlaringParameter[MAX_SPHERES],
    UniformVelocity[MAX_DIMENSION];
  FLOAT HaloRadius[MAX_SPHERES],
    HaloCoreRadius[MAX_SPHERES],
    HaloPosition[MAX_SPHERES][MAX_DIMENSION],
    DiskRadius[MAX_SPHERES],
    DiskHeight[MAX_SPHERES];
  float HaloMagneticField = 0.0;

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    HaloRadius[sphere]     = 1.0;
    HaloCoreRadius[sphere] = 0.1;
    HaloDensity[sphere]    = 1.0;
    HaloSpin[sphere]       = 0.05;
    HaloTemperature[sphere] = 1.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      HaloPosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] +
						   DomainRightEdge[dim]);
      HaloVelocity[sphere][dim] = 0;
    }
    GalaxyType[sphere]       = 0;
    DiskMassFraction[sphere] = 0.;
    DiskFlaringParameter[sphere] = 10.;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    UniformVelocity[dim] = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "NumberOfHalos = %"ISYM,
		  &NumberOfHalos);
    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "UseParticles = %"ISYM, 
		  &UseParticles);
    ret += sscanf(line, "UseGas = %"ISYM, 
		  &UseGas);
    ret += sscanf(line, "MediumTemperature = %"FSYM, 
		  &MediumTemperature);
    ret += sscanf(line, "MediumDensity = %"FSYM,
		  &MediumDensity);
    ret += sscanf(line, "HaloMagneticField = %"FSYM,
		  &HaloMagneticField);
    ret += sscanf(line, "UniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
		  UniformVelocity, UniformVelocity+1,
		  UniformVelocity+2);
    if (sscanf(line, "GalaxyType[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "GalaxyType[%"ISYM"] = %"ISYM, &sphere,
		    &GalaxyType[sphere]);
    if (sscanf(line, "HaloRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloRadius[%"ISYM"] = %"PSYM, &sphere,
		    &HaloRadius[sphere]);
    if (sscanf(line, "HaloCoreRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloCoreRadius[%"ISYM"] = %"PSYM, &sphere,
		    &HaloCoreRadius[sphere]);
    if (sscanf(line, "HaloDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloDensity[%"ISYM"] = %"FSYM, &sphere,
		    &HaloDensity[sphere]);
    if (sscanf(line, "HaloTemperature[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloTemperature[%"ISYM"] = %"FSYM, &sphere,
		    &HaloTemperature[sphere]);
    if (sscanf(line, "HaloAngVel[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloAngVel[%"ISYM"] = %"FSYM, &sphere,
		    &HaloAngVel[sphere]);
    if (sscanf(line, "HaloSpin[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloSpin[%"ISYM"] = %"FSYM, &sphere,
		    &HaloSpin[sphere]);
    if (sscanf(line, "HaloPosition[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloPosition[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &HaloPosition[sphere][0],
		    &HaloPosition[sphere][1],
		    &HaloPosition[sphere][2]);
    if (sscanf(line, "HaloVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "HaloVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &HaloVelocity[sphere][0],
		    &HaloVelocity[sphere][1],
		    &HaloVelocity[sphere][2]);
    if (sscanf(line, "DiskRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskRadius[%"ISYM"] = %"PSYM, &sphere,
		    &DiskRadius[sphere]);
    if (sscanf(line, "DiskHeight[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskHeight[%"ISYM"] = %"FSYM, &sphere,
		    &DiskHeight[sphere]);
    if (sscanf(line, "DiskDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskDensity[%"ISYM"] = %"FSYM, &sphere,
		    &DiskDensity[sphere]);
    if (sscanf(line, "DiskTemperature[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskTemperature[%"ISYM"] = %"FSYM, &sphere,
		    &DiskTemperature[sphere]);

    if (sscanf(line, "DiskMassFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskMassFraction[%"ISYM"] = %"FSYM, &sphere,
		    &DiskMassFraction[sphere]);

    if (sscanf(line, "DiskFlaringParameter[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "DiskFlaringParameter[%"ISYM"] = %"FSYM, &sphere,
		    &DiskFlaringParameter[sphere]);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxyDisk") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

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

  printf("timeu=%g(year)\n", TimeUnits/3.1558e7);
  printf("temp=%g, diskmassf=%g diskflarF=%g radius=%g, height=%g, density=%g\n",
	 DiskTemperature[0], DiskMassFraction[0], DiskFlaringParameter[0], DiskRadius[0], DiskHeight[0], DiskDensity[0]);

  if (UsePhysicalUnit) {
    MediumDensity /= DensityUnits;
    for (int sphere = 0; sphere < NumberOfHalos; sphere++) {
      HaloDensity[sphere] /= DensityUnits;
      HaloAngVel[sphere] *= TimeUnits;
      DiskDensity[sphere] /= DensityUnits;
    }
    HaloMagneticField /= MagneticUnits;
    printf("halodensity=%"GSYM"\n", HaloDensity[0]);
  }

  HaloVirialRadius  = HaloRadius[0]*LengthUnits;
  HaloConcentration = HaloVirialRadius/HaloCoreRadius[0]/LengthUnits;
  HaloCentralDensity = HaloDensity[0]*DensityUnits;

  if (DiskTemperature[0] > 0 && EOSSoundSpeed <= 0) {
    double tgamma = Gamma;
    if (EOSType == 3)
      tgamma = 1.;
    double c_s = sqrt(tgamma/Mu/mh * kboltz * DiskTemperature[0])/VelocityUnits;
    printf("EOSSoundSpeed was not set.\n");
    printf("Setting EOSSoundSpeed based on DiskTemperature[0]=%g K to %g (%g in code units)\n",
	   DiskTemperature[0], c_s*VelocityUnits, c_s);
    EOSSoundSpeed = c_s;
  }

  /* set up grid */

  if (TopGrid.GridData->GalaxyDiskInitializeGrid(
	     NumberOfHalos, HaloRadius,
	     HaloCoreRadius, HaloDensity,
	     HaloTemperature,
	     HaloPosition, HaloSpin,
	     HaloVelocity, HaloAngVel, HaloMagneticField,
	     DiskRadius, DiskHeight, 
	     DiskDensity, DiskTemperature, DiskMassFraction, DiskFlaringParameter, 
	     GalaxyType, UseParticles,
	     UseGas,
             UniformVelocity,
             MediumTemperature, MediumDensity, 0) == FAIL) {
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
	if (Temp->GridData->GalaxyDiskInitializeGrid(
	     NumberOfHalos, HaloRadius,
	     HaloCoreRadius, HaloDensity,
	     HaloTemperature,
	     HaloPosition, HaloSpin,
	     HaloVelocity, HaloAngVel, HaloMagneticField,
	     DiskRadius, DiskHeight, 
	     DiskDensity, DiskTemperature, DiskMassFraction, DiskFlaringParameter, 
	     GalaxyType, UseParticles,
	     UseGas,
	     UniformVelocity,
	     MediumTemperature, MediumDensity, level+1) == FAIL) {
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
  if(UsePoissonDivergenceCleaning){
    DataLabel[count++] = Phi_pName;
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

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  return SUCCESS;

}
