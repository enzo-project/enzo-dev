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

int Collapse1DInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";

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
  float SphereDensity,
    SpherePressure,
    SphereSoundVelocity,
    SphereAngVel;
  FLOAT SphereRadius,
    SphereCoreRadius;

  SphereRadius     = 1.0;
  SphereCoreRadius = 0.0;
  SphereDensity    = 1.0;
  SpherePressure   = 1.0;
  SphereSoundVelocity = 1.0;
  SphereAngVel = 0.0;
  SphereType       = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "RefineAtStart = %"ISYM, 
		  &RefineAtStart);
    ret += sscanf(line, "UseParticles = %"ISYM, 
		  &UseParticles);
    ret += sscanf(line, "MediumDensity = %"FSYM, 
		  &MediumDensity);
    ret += sscanf(line, "MediumPressure = %"FSYM,
		  &MediumPressure);

    ret += sscanf(line, "SphereType = %"ISYM,
		  &SphereType);
    ret += sscanf(line, "SphereRadius = %"PSYM,
		  &SphereRadius);
    ret += sscanf(line, "SphereCoreRadius = %"PSYM,
		  &SphereCoreRadius);
    ret += sscanf(line, "SphereDensity = %"FSYM,
		  &SphereDensity);
    ret += sscanf(line, "SpherePressure = %"FSYM,
		  &SpherePressure);
    ret += sscanf(line, "SphereSoundVelocity = %"FSYM,
		  &SphereSoundVelocity);
    ret += sscanf(line, "SphereAngVel = %"FSYM,
		  &SphereAngVel);
    /* if the line is suspicious, issue a warning */

  } // end input from parameter file
  
  //printf("InitialFractionHII=%"FSYM"\n", InitialFractionHII);
  //printf("1Frac = %"FSYM", Temp = %"FSYM"\n", FracKeplarianRot[sphere], SphereTemperature[sphere]);    

  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;
  }
  

  printf("rhou=%"GSYM",velu=%"GSYM",lenu=%"GSYM",tu=%"GSYM",presu=%"GSYM"\n", rhou, velu,lenu,tu,presu);

  // Bonnor-Ebert sphere: only the sound velocity and sphere radius are free parameters
  if (SphereType == 3) { 
    double G = 6.67e-8;
    double pi = 4.0*atan(1.0);

    double f=1.3; // BE sphere overdensity parameter
    double re = SphereRadius * lenu;
    double cs = SphereSoundVelocity;
    double ksi_e = 6.451; // critical radius of BE sphere
    double rhoc = ksi_e*ksi_e*f*cs*cs/(re*re*4*pi*G);

    SphereDensity = rhoc;
    MediumDensity = rhoc/14.0/100.0;
    MediumPressure = rhoc*cs*cs/14.0;
    double m_be = 1.18*pow(cs,4)/pow(G,1.5)/sqrt(MediumPressure);
    double msun = 1.989e33;
    m_be /= msun;

    printf("rhoc=%"GSYM", cs=%"GSYM", re=%"GSYM", m=%"GSYM"\n", rhoc, cs, re, m_be);
  }

  printf("rhoc=%"GSYM", rhom=%"GSYM", pm=%"GSYM"\n", SphereDensity, MediumDensity, MediumPressure);

  MediumDensity /= rhou;
  MediumPressure /= presu;

  SphereDensity /= rhou;
  SpherePressure /= presu;
  SphereSoundVelocity /= velu;
  SphereAngVel *= tu;

  if (TopGrid.GridData->Collapse1DInitializeGrid(
	     SphereRadius,
	     SphereCoreRadius, SphereDensity,
	     SpherePressure, SphereSoundVelocity, SphereAngVel, SphereType,
             MediumDensity, MediumPressure) == FAIL) {
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
      printf("In level %"ISYM"\n", level);
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->Collapse1DInitializeGrid(
				SphereRadius,
				SphereCoreRadius, SphereDensity,
				SpherePressure, SphereSoundVelocity, SphereAngVel, SphereType,
				MediumDensity, MediumPressure) == FAIL) {
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


  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  /* Write parameters to parameter output file */

  /*if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "NumberOfSpheres    = %"ISYM"\n",
	    n_sphere);
    fprintf(Outfptr, "RefineAtStart      = %"ISYM"\n",
	    RefineAtStart);
    fprintf(Outfptr, "UseParticles       = %"ISYM"\n",
	    UseParticles);
    fprintf(Outfptr, "UseColour          = %"ISYM"\n",
	    UseColour);
    fprintf(Outfptr, "InitialTemperature = %"FSYM"\n",
	    InitialTemperature);
    fprintf(Outfptr, "UniformVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
	    UniformVelocity[0], UniformVelocity[1],
	    UniformVelocity[2]);
    fprintf(Outfptr, "LengthUnit = %"FSYM"\n",
            LengthUnit);
    fprintf(Outfptr, "DensityUnit = %"FSYM"\n",
            DensityUnit);
    for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
      fprintf(Outfptr, "SphereType[%"ISYM"] = %"ISYM"\n", sphere,
	      SphereType[sphere]);
      fprintf(Outfptr, "SphereRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      SphereRadius[sphere]);
      fprintf(Outfptr, "SphereCoreRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      SphereCoreRadius[sphere]);
      fprintf(Outfptr, "SphereDensity[%"ISYM"] = %"FSYM"\n", sphere,
	      SphereDensity[sphere]);
      fprintf(Outfptr, "SphereTemperature[%"ISYM"] = %"FSYM"\n", sphere,
	      SphereTemperature[sphere]);
      fprintf(Outfptr, "SpherePosition[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			SpherePosition[sphere]);
      fprintf(Outfptr, "SphereVelocity[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			SphereVelocity[sphere]);
      fprintf(Outfptr, "FracKeplarianRot[%"ISYM"] = %"GOUTSYM"\n", sphere,
              FracKeplarianRot[sphere]);
      fprintf(Outfptr, "SphereTurbulence[%"ISYM"] = %"GOUTSYM"\n", sphere,
              SphereTurbulence[sphere]);
      fprintf(Outfptr, "SphereCutOff[%"ISYM"] = %"GOUTSYM"\n", sphere,
              SphereCutOff[sphere]);
      fprintf(Outfptr, "SphereAng1[%"ISYM"] = %"GOUTSYM"\n", sphere,
              SphereAng1[sphere]);
      fprintf(Outfptr, "SphereAng2[%"ISYM"] = %"GOUTSYM"\n", sphere,
              SphereAng2[sphere]);
      fprintf(Outfptr, "SphereNumShells[%"ISYM"] = %"ISYM"\n", sphere,
              SphereNumShells[sphere]);
    }
    }*/

  return SUCCESS;

}
