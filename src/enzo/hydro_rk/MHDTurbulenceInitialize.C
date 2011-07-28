/***********************************************************************
/
/  INITIALIZE MAGNETIZED TURBULENT CLOUD
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

int MHDTurbulenceInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData, int SetBaryonFields)
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
  char *Drive1Name = "DrivingField1";
  char *Drive2Name = "DrivingField2";
  char *Drive3Name = "DrivingField3";

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = TRUE;
  int RandomSeed = 1;
  float rho_medium=1.0, cs=1.0, mach=1.0, Bnaught=0.0;

  /* read input from file */
  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "RefineAtStart = %"ISYM, &RefineAtStart);
    ret += sscanf(line, "Density = %"FSYM, &rho_medium);
    ret += sscanf(line, "SoundVelocity = %"FSYM, &cs);
    ret += sscanf(line, "MachNumber = %"FSYM, &mach);
    ret += sscanf(line, "InitialBfield = %"FSYM, &Bnaught);
    ret += sscanf(line, "RandomSeed = %"ISYM, &RandomSeed);

  } // end input from parameter file

  /* Convert to code units */
  if (HydroMethod != MHD_RK) {// if its not with B-field make sure it is zero
    if (Bnaught != 0.) 
      printf("InitialBField was non-zero in paramterfile but hydromethod does not support a magnetic field. Proceeding with setting the initial field to zero.\n");
    Bnaught = 0.0;
  }

  printf(" RAW:  rho_medium = %"GSYM",cs = %"GSYM", mach = %"GSYM", Bnaught = %"GSYM" \n",rho_medium,cs,mach,Bnaught);

  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, 
    presu = 1.0, bfieldu = 1.0;
  if (UsePhysicalUnit) 
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
  presu = rhou*lenu*lenu/tu/tu;
  bfieldu = sqrt(presu*4.0*M_PI);
    
  rho_medium /= rhou;
  cs /= velu;
  Bnaught /= bfieldu;

  printf("rhou=%"GSYM",velu=%"GSYM",lenu=%"GSYM",tu=%"GSYM",presu=%"GSYM",bfieldu=%"GSYM", tempu=%"GSYM"\n", 
	 rhou, velu,lenu,tu,presu,bfieldu, tempu);
  printf("rho_medium=%"GSYM", cs=%"GSYM", Bnaught=%"GSYM"\n", rho_medium, cs, Bnaught);


  HierarchyEntry *CurrentGrid;

  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
    if (CurrentGrid->GridData->MHDTurbulenceInitializeGrid(rho_medium, cs, mach, 
				   Bnaught, RandomSeed, 0, SetBaryonFields) == FAIL) {
      fprintf(stderr, "Error in MHDTurbulenceInitializeGrid.\n");
      return FAIL;
    }

    CurrentGrid = CurrentGrid->NextGridThisLevel;
  }

  if (SetBaryonFields) {
    // Compute Normalization
    double v_rms  = 0;
    double Volume = 0;
    Eflt fac = 1;    

    CurrentGrid = &TopGrid;
    while (CurrentGrid != NULL) {
      if (CurrentGrid->GridData->PrepareVelocityNormalization(&v_rms, &Volume) == FAIL) {
	fprintf(stderr, "Error in PrepareVelocityNormalization.\n");
	return FAIL;
      }
      CurrentGrid = CurrentGrid->NextGridThisLevel;
      fprintf(stderr, "v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
    }
    
#ifdef USE_MPI
    CommunicationAllReduceValues(&v_rms, 1, MPI_SUM);
    CommunicationAllReduceValues(&Volume, 1, MPI_SUM);
#endif
    fprintf(stderr, "v_rms, Volume: %"GSYM"  %"GSYM"\n", v_rms, Volume);
    // Carry out the Normalization

    // Normalize Velocities now
    v_rms = sqrt(v_rms/Volume); // actuall v_rms
    fac = cs*mach/v_rms;

    CurrentGrid = &TopGrid;
    while (CurrentGrid != NULL) {
      if (CurrentGrid->GridData->NormalizeVelocities(fac) == FAIL) {
	fprintf(stderr, "Error in grid::NormalizeVelocities.\n");
	return FAIL;
      }
      CurrentGrid = CurrentGrid->NextGridThisLevel;
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
	if (Temp->GridData->MHDTurbulenceInitializeGrid(rho_medium, cs, fac, 
							Bnaught, RandomSeed, level, SetBaryonFields) == FAIL) {
	  fprintf(stderr, "Error in MHDTurbulenceInitializeGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels


    // Normalize Velocities now
    v_rms = sqrt(v_rms/Volume); // actuall v_rms
    fac = cs*mach/v_rms;

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->NormalizeVelocities(fac) == FAIL) {
	  fprintf(stderr, "Error in grid::NormalizeVelocities.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

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
  if (UseDrivingField) {
    DataLabel[count++] = Drive1Name;
    DataLabel[count++] = Drive2Name;
    DataLabel[count++] = Drive3Name;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  } // endif SetBaryonFields

  return SUCCESS;

}
