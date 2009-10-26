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

int CollapseMHD3DInitialize(FILE *fptr, FILE *Outfptr, 
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
  char *GravPotenName = "PotentialField";
  char *Acce1Name = "AccelerationField1";
  char *Acce2Name = "AccelerationField2";
  char *Acce3Name = "AccelerationField3";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int n_sphere = 1;
  int RefineAtStart   = TRUE;
  int UseParticles    = FALSE;
  float MediumDensity = 1.0, 
    MediumPressure = 1.0;
  int   SphereType[MAX_SPHERES];
  float SphereDensity[MAX_SPHERES],
    SpherePressure[MAX_SPHERES],
    SphereSoundVelocity[MAX_SPHERES],
    SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
    UniformVelocity[MAX_DIMENSION],
    SphereAngVel[MAX_SPHERES],
    SphereTurbulence[MAX_SPHERES],
    SphereCutOff[MAX_SPHERES],
    SphereAng1[MAX_SPHERES],
    SphereAng2[MAX_SPHERES];
  int	SphereNumShells[MAX_SPHERES];
  FLOAT SphereRadius[MAX_SPHERES],
    SphereCoreRadius[MAX_SPHERES],
    SpherePosition[MAX_SPHERES][MAX_DIMENSION];
  float Bnaught = 0.0;
  float theta_B = 0.0;

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    SphereRadius[sphere]     = 1.0;
    SphereCoreRadius[sphere] = 0.0;
    SphereDensity[sphere]    = 1.0;
    SpherePressure[sphere]   = 1.0;
    SphereSoundVelocity[sphere] = 1.0;
    SphereAngVel[sphere] = 0.0;
    SphereTurbulence[sphere] = 0.0;
    SphereCutOff[sphere] = 6.5;
    SphereAng1[sphere] = 0;
    SphereAng2[sphere] = 0;
    SphereNumShells[sphere] = 1;

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SpherePosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
      SphereVelocity[sphere][dim] = 0;
    }
    SphereType[sphere]       = 0;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    UniformVelocity[dim] = 0;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "NumberOfSpheres = %d", &n_sphere);
    ret += sscanf(line, "RefineAtStart = %d", &RefineAtStart);
    ret += sscanf(line, "UseParticles = %d", &UseParticles);
    ret += sscanf(line, "MediumDensity = %f", &MediumDensity);
    ret += sscanf(line, "MediumPressure = %f", &MediumPressure);
    ret += sscanf(line, "UniformVelocity = %f %f %f", 
		  UniformVelocity, UniformVelocity+1,
		  UniformVelocity+2);
    ret += sscanf(line, "InitialBField = %f", &Bnaught);
    ret += sscanf(line, "theta_B = %f", &theta_B);
 
    if (sscanf(line, "SphereType[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereType[%d] = %d", &sphere,
		    &SphereType[sphere]);
    if (sscanf(line, "SphereRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereRadius[%d] = %"FSYM, &sphere,
		    &SphereRadius[sphere]);
    if (sscanf(line, "SphereCoreRadius[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereCoreRadius[%d] = %"FSYM, &sphere,
		    &SphereCoreRadius[sphere]);
    if (sscanf(line, "SphereDensity[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereDensity[%d] = %f", &sphere,
		    &SphereDensity[sphere]);
    if (sscanf(line, "SpherePressure[%d]", &sphere) > 0)
      ret += sscanf(line, "SpherePressure[%d] = %f", &sphere,
		    &SpherePressure[sphere]);
    if (sscanf(line, "SphereSoundVelocity[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereSoundVelocity[%d] = %f", &sphere,
		    &SphereSoundVelocity[sphere]);
    if (sscanf(line, "SpherePosition[%d]", &sphere) > 0)
      ret += sscanf(line, "SpherePosition[%d] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &SpherePosition[sphere][0],
		    &SpherePosition[sphere][1],
		    &SpherePosition[sphere][2]);
    if (sscanf(line, "SphereVelocity[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereVelocity[%d] = %f %f %f", 
		    &sphere, &SphereVelocity[sphere][0],
		    &SphereVelocity[sphere][1],
		    &SphereVelocity[sphere][2]);
    if (sscanf(line, "SphereAngVel[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereAngVel[%d] = %f", &sphere,
                    &SphereAngVel[sphere]);
    if (sscanf(line, "SphereTurbulence[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereTurbulence[%d] = %f", &sphere,
                    &SphereTurbulence[sphere]);
    if (sscanf(line, "SphereCutOff[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereCutOff[%d] = %f", &sphere,
                    &SphereCutOff[sphere]);
    if (sscanf(line, "SphereAng1[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereAng1[%d] = %f", &sphere,
                    &SphereAng1[sphere]);
    if (sscanf(line, "SphereAng2[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereAng2[%d] = %f", &sphere,
                    &SphereAng2[sphere]);
    if (sscanf(line, "SphereNumShells[%d]", &sphere) > 0)
      ret += sscanf(line, "SphereNumShells[%d] = %d", &sphere,
                    &SphereNumShells[sphere]);
    /* if the line is suspicious, issue a warning */

  } // end input from parameter file
  
  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0, bfieldu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;
    bfieldu = sqrt(presu*4.0*M_PI);
  }
  
  printf("rhou=%g,velu=%g,lenu=%g,tu=%g (%g yr),tempu=%g,presu=%g, bfieldu=%g, tempu=%g\n", 
	 rhou, velu,lenu,tu,tu/3.1558e7,tempu,presu,bfieldu, tempu);

  // Bonnor-Ebert sphere: only the sound velocity and sphere radius are free parameters
  if (SphereType[0] == 3) { 
    double G = 6.67e-8;

    double f=1.5; // BE sphere overdensity parameter
    double re = SphereRadius[0] * lenu;
    double cs = SphereSoundVelocity[0];
    double ksi_e = 6.451; // critical radius of BE sphere
    double rhoc = ksi_e*ksi_e*f*cs*cs/(re*re*4*M_PI*G);

    SphereDensity[0] = rhoc;
    //    MediumDensity = rhoc/14.0;
    
    MediumPressure = rhoc*cs*cs/14.0;
    double m_be = pow(f,1.5)*1.18*pow(cs,4)/pow(G,1.5)/sqrt(MediumPressure);
    double msun = 1.989e33;
    m_be /= msun;

    printf("rhoc=%g, cs=%g, re=%g, m=%g\n", rhoc, cs, re, m_be);
  }


  MediumDensity /= rhou;
  MediumPressure /= presu;

  Bnaught /= bfieldu;

  //printf("t=%g\n", MediumPressure/MediumDensity*tempu);

  for (int i = 0; i < n_sphere; i++) {
    SphereDensity[i] /= rhou;
    SpherePressure[i] /= presu;
    SphereSoundVelocity[i] /= velu;
    SphereAngVel[i] *= tu;
  }

  printf("rhoc=%g, rhom=%g, pm=%g\n", SphereDensity[0], MediumDensity, MediumPressure);


  if (TopGrid.GridData->CollapseMHD3DInitializeGrid(
	     n_sphere, SphereRadius,
	     SphereCoreRadius, SphereDensity,
	     SpherePressure, SphereSoundVelocity, SpherePosition, SphereAngVel, Bnaught, theta_B,
	     SphereType,
             MediumDensity, MediumPressure, 0) == FAIL) {
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
				n_sphere, SphereRadius,
				SphereCoreRadius, SphereDensity,
				SpherePressure, SphereSoundVelocity, SpherePosition, SphereAngVel, Bnaught,theta_B,
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

  if (WritePotential) {
    DataLabel[count++] = GravPotenName;
    DataLabel[count++] = Acce1Name;
    DataLabel[count++] = Acce2Name;
    DataLabel[count++] = Acce3Name;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return SUCCESS;

}
