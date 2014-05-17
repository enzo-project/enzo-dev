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
#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

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
#include "CommunicationUtilities.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CollapseMHD3DInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData, int SetBaryonFields)
{
   const char *DensName = "Density";
   const char *TEName   = "TotalEnergy";
   const char *GEName   = "GasEnergy";
   const char *Vel1Name = "x-velocity";
   const char *Vel2Name = "y-velocity";
   const char *Vel3Name = "z-velocity";
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
   const char *DebugName = "Debug";
   const char *Phi_pName = "Phip";
   const char *GravPotenName = "PotentialField";
   const char *Acce1Name = "AccelerationField1";
   const char *Acce2Name = "AccelerationField2";
   const char *Acce3Name = "AccelerationField3";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int n_sphere = 1;
  int RefineAtStart   = TRUE;
  int UseParticles    = FALSE;
  float MediumDensity = 1.0, 
    MediumPressure = 1.0;
  float MediumTemperature = FLOAT_UNDEFINED;
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
  int Bdirection = 2;

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

    ret += sscanf(line, "NumberOfSpheres = %"ISYM, &n_sphere);
    ret += sscanf(line, "RefineAtStart = %"ISYM, &RefineAtStart);
    ret += sscanf(line, "UseParticles = %"ISYM, &UseParticles);
    ret += sscanf(line, "MediumDensity = %"FSYM, &MediumDensity);
    ret += sscanf(line, "MediumPressure = %"FSYM, &MediumPressure);
    ret += sscanf(line, "MediumTemperature = %"FSYM, &MediumTemperature);
    ret += sscanf(line, "UniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
		  UniformVelocity, UniformVelocity+1,
		  UniformVelocity+2);
    ret += sscanf(line, "InitialBField = %"FSYM, &Bnaught);
    ret += sscanf(line, "theta_B = %"FSYM, &theta_B);
    ret += sscanf(line, "Bdirection = %"ISYM, &Bdirection);

    if (sscanf(line, "SphereType[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereType[%"ISYM"] = %"ISYM, &sphere,
		    &SphereType[sphere]);
    if (sscanf(line, "SphereRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereRadius[%"ISYM"] = %"PSYM, &sphere,
		    &SphereRadius[sphere]);
    if (sscanf(line, "SphereCoreRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereCoreRadius[%"ISYM"] = %"PSYM, &sphere,
		    &SphereCoreRadius[sphere]);
    if (sscanf(line, "SphereDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereDensity[%"ISYM"] = %"FSYM, &sphere,
		    &SphereDensity[sphere]);
    if (sscanf(line, "SpherePressure[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SpherePressure[%"ISYM"] = %"FSYM, &sphere,
		    &SpherePressure[sphere]);
    if (sscanf(line, "SphereSoundVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereSoundVelocity[%"ISYM"] = %"FSYM, &sphere,
		    &SphereSoundVelocity[sphere]);
    if (sscanf(line, "SpherePosition[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SpherePosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM, 
		    &sphere, &SpherePosition[sphere][0],
		    &SpherePosition[sphere][1],
		    &SpherePosition[sphere][2]);
    if (sscanf(line, "SphereVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &SphereVelocity[sphere][0],
		    &SphereVelocity[sphere][1],
		    &SphereVelocity[sphere][2]);
    if (sscanf(line, "SphereAngVel[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereAngVel[%"ISYM"] = %"FSYM, &sphere,
                    &SphereAngVel[sphere]);
    if (sscanf(line, "SphereTurbulence[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereTurbulence[%"ISYM"] = %"FSYM, &sphere,
                    &SphereTurbulence[sphere]);
    if (sscanf(line, "SphereCutOff[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereCutOff[%"ISYM"] = %"FSYM, &sphere,
                    &SphereCutOff[sphere]);
    if (sscanf(line, "SphereAng1[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereAng1[%"ISYM"] = %"FSYM, &sphere,
                    &SphereAng1[sphere]);
    if (sscanf(line, "SphereAng2[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereAng2[%"ISYM"] = %"FSYM, &sphere,
                    &SphereAng2[sphere]);
    if (sscanf(line, "SphereNumShells[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "SphereNumShells[%"ISYM"] = %"ISYM, &sphere,
                    &SphereNumShells[sphere]);
    /* if the line is suspicious, issue a warning */

  } // end input from parameter file
  
  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0, bfieldu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;
    bfieldu = sqrt(presu*4.0*M_PI);
  }
  
  printf("rhou=%"GSYM",velu=%"GSYM",lenu=%"GSYM",tu=%"GSYM" (%"GSYM" yr),tempu=%"GSYM",presu=%"GSYM", bfieldu=%"GSYM", tempu=%"GSYM"\n", 
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
    MediumDensity *= rhoc/14.0; // in this case medium density is the ratio to decrease the outside density 
    
    MediumPressure = rhoc/14.0 *cs*cs/Gamma;
    double m_be = pow(f,1.5)*1.18*pow(cs,4)/pow(G,1.5)/sqrt(MediumPressure);
    double msun = 1.989e33;
    m_be /= msun;

    printf("rhoc=%"GSYM", cs=%"GSYM", re=%"GSYM", m=%"GSYM"\n", rhoc, cs, re, m_be);
  }


  MediumDensity /= rhou;
  if (MediumTemperature > 0) {
    float MediumEnergy;
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
      VelocityUnits;
    float tmp1, tmp2, tmp3;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	     &VelocityUnits, MetaData.Time);
    MediumEnergy = MediumTemperature / TemperatureUnits / ((Gamma-1.0)*Mu);
    MediumPressure = (Gamma-1.0) * MediumDensity * MediumEnergy;
  } else {
    MediumPressure /= presu;
  }

  Bnaught /= bfieldu;

  //printf("t=%"GSYM"\n", MediumPressure/MediumDensity*tempu);

  for (int i = 0; i < n_sphere; i++) {
    SphereDensity[i] /= rhou;
    SpherePressure[i] /= presu;
    SphereSoundVelocity[i] /= velu;
    SphereAngVel[i] *= tu;
  }

  printf("rhoc=%"GSYM", rhom=%"GSYM", pm=%"GSYM"\n", SphereDensity[0], MediumDensity, MediumPressure);

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  int count = 0;
  while (CurrentGrid != NULL) {
    printf("count %i %i\n", count++, MyProcessorNumber);
    if (CurrentGrid->GridData->CollapseMHD3DInitializeGrid(
						      n_sphere, SphereRadius,
						      SphereCoreRadius, SphereDensity,
						      SpherePressure, SphereSoundVelocity, SpherePosition, 
						      SphereAngVel, SphereTurbulence, Bnaught, theta_B, Bdirection,
						      SphereType,
						      MediumDensity, MediumPressure, 0, SetBaryonFields) == FAIL) {
      fprintf(stderr, "Error in CollapseMHD3DInitializeGrid.\n");
      return FAIL;
    }
    CurrentGrid = CurrentGrid->NextGridThisLevel;
  }


  if (SetBaryonFields) {

    // Compute Velocity Normalization
    double v_rms  = 0;
    double Volume = 0;
    Eflt fac = 1;
    
    if (SphereTurbulence[0] > 0.) {
      CurrentGrid = &TopGrid;
      while (CurrentGrid != NULL) {
	if (CurrentGrid->GridData->PrepareVelocityNormalization(&v_rms, &Volume) == FAIL) {
	  fprintf(stderr, "Error in PrepareVelocityNormalization.\n");
	  return FAIL;
	}
	CurrentGrid = CurrentGrid->NextGridThisLevel;
	fprintf(stderr, "Prepared: v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
      }
      
#ifdef USE_MPI
      CommunicationAllReduceValues(&v_rms, 1, MPI_SUM);
      CommunicationAllReduceValues(&Volume, 1, MPI_SUM);
#endif
      fprintf(stderr, "v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
      // Carry out the Normalization
      v_rms = sqrt(v_rms/Volume); // actuall v_rms
      fac = SphereSoundVelocity[0]*SphereTurbulence[0]/v_rms;
      
      CurrentGrid = &TopGrid;
      while (CurrentGrid != NULL) {
	if (CurrentGrid->GridData->NormalizeVelocities(fac) == FAIL) {
	  fprintf(stderr, "Error in grid::NormalizeVelocities.\n");
	  return FAIL;
	}
	CurrentGrid = CurrentGrid->NextGridThisLevel;
      }
      if (fac != 0. ) SphereTurbulence[0] = fac;
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
							  SpherePressure, SphereSoundVelocity, 
							  SpherePosition, SphereAngVel, SphereTurbulence,  
							  Bnaught, theta_B, Bdirection,
							  SphereType, MediumDensity, MediumPressure, level+1,
							  SetBaryonFields) == FAIL) {
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

  } // end if  SetBaryonField

    /* set up field names and units */
  count = 0;
  DataLabel[count++] = (char*)  DensName;
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
  if (UseMHDCT) {
      MHDcLabel[0] = "Bx";
      MHDcLabel[1] = "By";
      MHDcLabel[2] = "Bz";

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
  
  if(UseDivergenceCleaning){
    DataLabel[count++] = (char*) Phi_pName;
    DataLabel[count++] = (char*) DebugName;
  }
  
  if (WritePotential) {
    DataLabel[count++] = (char*) GravPotenName;
    //    DataLabel[count++] = (char*) Acce1Name;
    //    DataLabel[count++] = (char*) Acce2Name;
    //    DataLabel[count++] = (char*) Acce3Name;
  }
  
  DataLabel[count++] = (char*) DebugName;
  
  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }
  
  
  return SUCCESS;

}
