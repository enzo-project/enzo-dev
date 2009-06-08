/***********************************************************************
/
/  INITIALIZE COOLING TEST
/
/  written by: John Wise
/  date:       April 2009
/  modified1:
/
/  PURPOSE:
/    Set up a grid that has density, temperature, and a species 
/    (e.g. metallicity) log-spaced in each direction.  Original idea 
/    taken from Britton Smith.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//
#include <stdlib.h>
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

int CoolingTestInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
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
  const char *MetalName = "Metal_Density";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i, source;

  /* set default parameters */

  float CoolingTestMinimumDensity = 1e-2;
  float CoolingTestMaximumDensity = 1e6;
  float CoolingTestMinimumTemperature = 10.0;
  float CoolingTestMaximumTemperature = 1e8;
  float CoolingTestMinimumColour = 1e-6;
  float CoolingTestMaximumColour = 1.0;
  float CoolingTestElectronFraction = 1e-3;
  float CoolingTestH2Fraction = 1e-3;
  int CoolingTestUseMetals = TRUE;
  int CoolingTestUseElectronFraction = FALSE;

  /* read input from file */
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CoolingTestMinimumDensity = %"FSYM,
		  &CoolingTestMinimumDensity);
    ret += sscanf(line, "CoolingTestMaximumDensity = %"FSYM,
		  &CoolingTestMaximumDensity);
    ret += sscanf(line, "CoolingTestMinimumTemperature = %"FSYM,
		  &CoolingTestMinimumTemperature);
    ret += sscanf(line, "CoolingTestMaximumTemperature = %"FSYM,
		  &CoolingTestMaximumTemperature);
    ret += sscanf(line, "CoolingTestMinimumColour = %"FSYM,
		  &CoolingTestMinimumColour);
    ret += sscanf(line, "CoolingTestMaximumColour = %"FSYM,
		  &CoolingTestMaximumColour);
    ret += sscanf(line, "CoolingTestElectronFraction = %"FSYM,
		  &CoolingTestElectronFraction);
    ret += sscanf(line, "CoolingTestH2Fraction = %"FSYM,
		  &CoolingTestH2Fraction);
    ret += sscanf(line, "CoolingTestUseMetals = %"ISYM,
		  &CoolingTestUseMetals);
    ret += sscanf(line, "CoolingTestUseElectronFraction = %"ISYM,
		  &CoolingTestUseElectronFraction);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CoolingTest") && 
	line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning62: %"ISYM", the following parameter line was "
		"not interpreted:\n%s\n", ret, line);
    
  } // end input from parameter file

  /* Error check */

  if (CoolingTestUseMetals && CoolingTestUseElectronFraction) {
    fprintf(stderr, 
	    "WARNING: CoolingTestUseMetals and CoolingTestUseElectronFraction "
	    "are both TRUE.  Please pick one or the other.\n");
    return FAIL;
  }

  /* Override some parameters and defaults */

  MetaData.StaticHierarchy = TRUE;
  HydroMethod = HydroMethodUndefined;
  if (CoolingTestUseElectronFraction)
    if (MultiSpecies < 1) {
      fprintf(stderr, 
	      "WARNING: MultiSpecies must be ON for CoolingTestUseElectronFraction.\n"
	      "         Turning ON.\n");
      MultiSpecies = 1;
    }

  if (CoolingTestUseMetals) {
    MetalCooling = JHW_METAL_COOLING;
    if (MultiSpecies < 1) {
      fprintf(stderr, 
	      "WARNING: MultiSpecies must be ON for CoolingTestUseMetals.\n"
	      "         Turning ON.\n");
      MultiSpecies = 1;
    }

    // Check for metal cooling rate table

    FILE *test_fptr = fopen(MetalCoolingTable, "r");
    if (test_fptr == NULL) {
      fprintf(stderr, "Error opening metal cooling table %s\n", MetalCoolingTable);
      return FAIL;
    }
    fclose(test_fptr);

    // Estimate species fractions from electron fraction

    TestProblemData.HI_Fraction = CoolData.HydrogenFractionByMass * 
      (1.0 - CoolingTestElectronFraction);
    TestProblemData.HII_Fraction = CoolData.HydrogenFractionByMass * 
      CoolingTestElectronFraction;
    TestProblemData.HeI_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) * 
      0.5*(1.0 - CoolingTestElectronFraction);
    TestProblemData.HeII_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) * 
      0.5*(1.0 - CoolingTestElectronFraction);
    TestProblemData.HeIII_Fraction = 4*(1.0 - CoolData.HydrogenFractionByMass) *
      CoolingTestElectronFraction;

  }

  TestProblemData.H2I_Fraction = CoolingTestH2Fraction;
  TestProblemData.H2II_Fraction = tiny_number;

  /* set up grid */

  if (TopGrid.GridData->
      CoolingTestInitializeGrid(CoolingTestMinimumDensity,
				CoolingTestMaximumDensity,
				CoolingTestMinimumTemperature,
				CoolingTestMaximumTemperature,
				CoolingTestMinimumColour,
				CoolingTestMaximumColour,
				CoolingTestUseMetals,
				CoolingTestUseElectronFraction) == FAIL) {
    fprintf(stderr, "Error in CoolingTestInitializeGrid.\n");
    return FAIL;
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = (char*) GEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
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
  if (CoolingTestUseMetals)
    DataLabel[count++] = (char*) MetalName;
  
  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CoolingTestMinimumDensity = %"GOUTSYM"\n",
		  CoolingTestMinimumDensity);
    fprintf(Outfptr, "CoolingTestMaximumDensity = %"GOUTSYM"\n",
		  CoolingTestMaximumDensity);
    fprintf(Outfptr, "CoolingTestMinimumTemperature = %"GOUTSYM"\n",
		  CoolingTestMinimumTemperature);
    fprintf(Outfptr, "CoolingTestMaximumTemperature = %"GOUTSYM"\n",
		  CoolingTestMaximumTemperature);
    fprintf(Outfptr, "CoolingTestMinimumColour = %"GOUTSYM"\n",
		  CoolingTestMinimumColour);
    fprintf(Outfptr, "CoolingTestMaximumColour = %"GOUTSYM"\n",
		  CoolingTestMaximumColour);
    fprintf(Outfptr, "CoolingTestUseMetals = %"ISYM"\n",
		  CoolingTestUseMetals);
    fprintf(Outfptr, "CoolingTestUseElectronFraction = %"ISYM"\n\n",
		  CoolingTestUseElectronFraction);
  }

  return SUCCESS;

}
