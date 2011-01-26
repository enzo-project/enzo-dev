/***********************************************************************
/
/  INITIALIZE A SUPERNOVA EXPLOSION FROM A RESTART CALCULATION
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:
/
/  PURPOSE:  This routine reads in a previously generated output file
/            (presumably from a cosmology calculation), and then
/            initializes a supernova explosion in the center.
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
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
 
int SupernovaRestartInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior)
{
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, level, ret;
  float dummyf;
  /* Set default supernova parameters. */
 
  float SupernovaRestartEjectaMass   = 1.0;   // in solar masses
  float SupernovaRestartEjectaRadius = 1.0;   // in pc
  float SupernovaRestartEjectaEnergy = 1.0;   // in 10^51 erg
  FLOAT SupernovaRestartEjectaCenter[MAX_DIMENSION];
  int   SupernovaRestartColourField   = FALSE;
 
  char *SupernovaRestartName = NULL;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    SupernovaRestartEjectaCenter[dim] = FLOAT_UNDEFINED;
 
  /* Error check. */
 
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* Read parameters */
 
    ret += sscanf(line, "SupernovaRestartEjectaMass = %"FSYM,
		  &SupernovaRestartEjectaMass);
    ret += sscanf(line, "SupernovaRestartEjectaRadius = %"FSYM,
		  &SupernovaRestartEjectaRadius);
    ret += sscanf(line, "SupernovaRestartEjectaEnergy = %"FSYM,
		  &SupernovaRestartEjectaEnergy);
    ret += sscanf(line, "SupernovaRestartEjectaCenter = %"PSYM" %"PSYM" %"PSYM,
		  SupernovaRestartEjectaCenter,
		  SupernovaRestartEjectaCenter+1,
		  SupernovaRestartEjectaCenter+2);
    ret += sscanf(line, "SupernovaRestartColourField = %"ISYM,
		  &SupernovaRestartColourField);
 
    if (sscanf(line, "SupernovaRestartName = %s", dummy) == 1)
      SupernovaRestartName = dummy;
 
    /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "SupernovaRestart") &&
	line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* More error checking. */
 
  if (SupernovaRestartName == NULL) {
    ENZO_FAIL("Missing restart file name.\n");
  }
 
  /* -------------------------------------------------------------------- */
  /* Read the restart file. */
 
  if (debug)
    printf("reading restart parameter file %s\n", SupernovaRestartName);
  if (ReadAllData(SupernovaRestartName, &TopGrid, MetaData, &Exterior, &dummyf)
      == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
    ENZO_VFAIL("Error in ParameterFile %s.\n", SupernovaRestartName)
  }
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(stderr, "Successfully read restart file %s.\n",
	    SupernovaRestartName);
 
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
 
  float EjectaRadius = SupernovaRestartEjectaRadius * LengthConversion;
  float EjectaDensity = SupernovaRestartEjectaMass * MassConversion/
                        (4.0/3.0*3.14159*POW(EjectaRadius, 3));
  float EjectaThermalEnergy = SupernovaRestartEjectaEnergy * EnergyConversion /
        (SupernovaRestartEjectaMass * MassConversion);
 
  EjectaRadius        /= LengthUnits;
  EjectaDensity       /= DensityUnits;
  EjectaThermalEnergy /= VelocityUnits*VelocityUnits;
 
  if (debug) {
    printf("SupernovaRestart: initial T = %"GSYM" K\n",
	   EjectaThermalEnergy*TemperatureUnits*(Gamma-1.0)*0.6);
    printf("SupernovaRestart: r (code units) = %"GSYM"\n", EjectaRadius);
    printf("SupernovaRestart: density (code units) = %"GSYM"\n", EjectaDensity);
  }
 
  /* -------------------------------------------------------------------- */
  /* Loop over all the grids and call the initializer to modify them
     if necessary. */
 
  LevelHierarchyEntry *Temp, *Temp2;
  int NumberOfCellsSet = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    Temp = LevelArray[level];
 
    while (Temp != NULL) {
 
      if (Temp->GridData->SupernovaRestartInitialize(EjectaDensity,
			       EjectaRadius, EjectaThermalEnergy,
			       SupernovaRestartEjectaCenter,
			       SupernovaRestartColourField,
			       &NumberOfCellsSet) == FAIL) {
	ENZO_FAIL("Error in grid->SupernovaRestartInitialize\n");
      }
      Temp = Temp->NextGridThisLevel;
    }
 
  }
  if (debug)
    printf("SupernovaRestart: NumberOfCellsSet = %"ISYM"\n", NumberOfCellsSet);
 
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

    fprintf(Outfptr, "SupernovaRestartEjectaMass   = %"FSYM"\n",
	    SupernovaRestartEjectaMass);
    fprintf(Outfptr, "SupernovaRestartEjectaRadius = %"FSYM"\n",
	    SupernovaRestartEjectaRadius);
    fprintf(Outfptr, "SupernovaRestartEjectaEnergy = %"FSYM"\n",
	    SupernovaRestartEjectaEnergy);
    fprintf(Outfptr, "SupernovaRestartEjectaCenter = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank,
		      SupernovaRestartEjectaCenter);
    fprintf(Outfptr, "SupernovaRestartColourField  = %"ISYM"\n",
	    SupernovaRestartColourField);
    fprintf(Outfptr, "SupernovaRestartName         = %s\n",
	    SupernovaRestartName);
 
  }
 
  /* Clean up. */
 
  delete dummy;
 
  return SUCCESS;
}
