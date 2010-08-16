/***********************************************************************
/
/  INITIALIZE PHOTON TEST FROM A RESTART CALCULATION
/
/  written by: Elizabeth Harper-Clark
/  date:       March, 2010
/  modified1:
/
/  PURPOSE:  Based on SuperNovaeRestartInitialize files
/            This routine reads in a previously generated output file
/            (presumably from a cosmology calculation), and then puts 
/            in Sink particles of specified mass and position.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
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
#include "CosmologyParameters.h"
#include "fortran.def"
 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		ExternalBoundary *Exterior, float *Initialdt);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		      ExternalBoundary *Exterior, float *Initialdt,
		      bool ReadParticlesOnly=false);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
 
int PhotonTestRestartInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior)
{
  printf("       PhotonTestRestart STARTED\n"); 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, level, ret;
  float dummyf;
 
  char *PhotonTestRestartName = NULL;
 
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    if (sscanf(line, "PhotonTestRestartName = %s", dummy) == 1)
      PhotonTestRestartName = dummy;
 
    /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "PhotonTestRestart") &&
	line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* More error checking. */
 
  if (PhotonTestRestartName == NULL) {
    ENZO_FAIL("Missing restart file name.\n");
  }

  /* -------------------------------------------------------------------- */
  /* Read the restart file. */
 
  if (debug)
    printf("reading restart parameter file %s\n", PhotonTestRestartName);
 
#ifdef USE_HDF5_GROUPS
  if (Group_ReadAllData(PhotonTestRestartName, &TopGrid, MetaData, &Exterior, &dummyf) == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Error in Group_ReadAllData %s\n",PhotonTestRestartName );
      fprintf(stderr, "Probably not in a packed-HDF5 format. Trying other read routines.\n");
    }
#endif
    // If not packed-HDF5, then try usual HDF5 or HDF4
    if (ReadAllData(PhotonTestRestartName, &TopGrid, MetaData, &Exterior, &dummyf)
	== FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
      ENZO_VFAIL("Error in ParameterFile %s.\n", PhotonTestRestartName)
    }
#ifdef USE_HDF5_GROUPS
  }
#endif

  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Successfully read restart file %s.\n",
	    PhotonTestRestartName);

  if (ProblemType == 51)
    if (ReadPhotonSources(fptr, MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in ReadPhotonSources.\n");
    }
 
  PhotonTime = InitialTimeInCodeUnits;

  printf("    PhotonTestRestartName = %s   \n",PhotonTestRestartName ); 

  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
 
  /* Convert Mass, Radius and Energy into code units of density, radius
     and thermal energy (by assuming mass and energy is evenly distributed
     within radius).
     If comoving coordinate is set then assume mass is in solar units,
     radius in pc, and energy in 10^51 erg; otherwise do no conversion.*/
 
 
  double MassConversion = 1, LengthConversion = 1, EnergyConversion = 1;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
    TemperatureUnits = 1;
 
  FLOAT Time = TopGrid.GridData->ReturnTime();
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  if (ComovingCoordinates) {
 
    LengthConversion = 3.08e18;     // pc
    MassConversion   = 2e33;        // solar masses
    EnergyConversion = 1.0e51;      // 10^51 erg
 
  }
 
//   float EjectaRadius = SupernovaRestartEjectaRadius * LengthConversion;
//   float EjectaDensity = SupernovaRestartEjectaMass * MassConversion/
//                         (4.0/3.0*3.14159*POW(EjectaRadius, 3));
//   float EjectaThermalEnergy = SupernovaRestartEjectaEnergy * EnergyConversion /
//         (SupernovaRestartEjectaMass * MassConversion);
 
//   EjectaRadius        /= LengthUnits;
//   EjectaDensity       /= DensityUnits;
//   EjectaThermalEnergy /= VelocityUnits*VelocityUnits;
 
//   if (debug) {
//     printf("PhotonTestRestart: initial T = %"GSYM" K\n",
// 	   EjectaThermalEnergy*TemperatureUnits*(Gamma-1.0)*0.6);
//     printf("PhotonTestRestart: r (code units) = %"GSYM"\n", EjectaRadius);
//     printf("PhotonTestRestart: density (code units) = %"GSYM"\n", EjectaDensity);
//   }
 
  /* -------------------------------------------------------------------- */
  /* Loop over all the grids and call the initializer to modify them
     if necessary. */
 
  LevelHierarchyEntry *Temp, *Temp2;
  int NumberOfCellsSet = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    Temp = LevelArray[level];
 
    while (Temp != NULL) {
 
      if (Temp->GridData->PhotonTestRestartInitialize(level,&NumberOfCellsSet) == FAIL) {
	ENZO_FAIL("Error in grid->PhotonTestRestartInitialize\n");
      }
      Temp = Temp->NextGridThisLevel;
    }
 
  }
  if (debug)
    printf("PhotonTestRestart: NumberOfCellsSet = %"ISYM"\n", NumberOfCellsSet);
 
  /* -------------------------------------------------------------------- */
  /* Loop over grid and project solution to parent to maintain consistency. */
 
  for (level = MaximumRefinementLevel; level > 0; level--) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridData->ProjectSolutionToParentGrid(
                *(Temp->GridHierarchyEntry->ParentGrid->GridData)) == FAIL) {
	ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");
      }
      Temp2 = Temp->NextGridThisLevel;
      delete Temp;   // clean up as we go along
      Temp = Temp2;
    }
  }
 
  /* -------------------------------------------------------------------- */
  /* Write parameters to parameter output file */
 
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "PhotonTestRestartName         = %s\n",
	    PhotonTestRestartName);
 
  }
 
  /* Clean up. */
 
  delete dummy;
 
  return SUCCESS;
}
