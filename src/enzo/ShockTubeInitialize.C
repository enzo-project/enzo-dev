/***********************************************************************
/
/  INITIALIZE A SHOCK TUBE SIMULATION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
 
int ShockTubeInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int ret;
  float ShockTubeDensity[2], ShockTubePressure[2], ShockTubeVelocity[2];
 
  /* set default parameters */
 
  int ShockTubeDirection = 0;      // x direction
  float ShockTubeBoundary = 0.5;     //
 
  ShockTubeDensity[0]  = 1.0;  // This is the classic Sod Shock tube problem
  ShockTubePressure[0] = 1.0;
  ShockTubeVelocity[0] = 0.0;
 
  ShockTubeDensity[1]  = 0.125;
  ShockTubePressure[1] = 0.1;
  ShockTubeVelocity[1] = 0.0;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "ShockTubeDirection = %"ISYM, &ShockTubeDirection);
    ret += sscanf(line, "ShockTubeBoundary = %"PSYM, &ShockTubeBoundary);
 
    ret += sscanf(line, "ShockTubeLeftDensity = %"PSYM, &ShockTubeDensity[0]);
    ret += sscanf(line, "ShockTubeLeftPressure = %"PSYM, &ShockTubePressure[0]);
    ret += sscanf(line, "ShockTubeLeftVelocity = %"PSYM, &ShockTubeVelocity[0]);
 
    ret += sscanf(line, "ShockTubeRightDensity = %"PSYM, &ShockTubeDensity[1]);
    ret += sscanf(line, "ShockTubeRightPressure = %"PSYM, &ShockTubePressure[1]);
    ret += sscanf(line, "ShockTubeRightVelocity = %"PSYM, &ShockTubeVelocity[1]);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "ShockTube"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */
 
  if (TopGrid.GridData->ShockTubeInitializeGrid(ShockTubeDirection,
						ShockTubeBoundary,
						ShockTubeDensity,
						ShockTubePressure,
						ShockTubeVelocity) == FAIL) {
    fprintf(stderr, "Error in ShockTubeInitializeGrid.\n");
    return FAIL;
  }
 
  /* set up field names and units */
 
  DataLabel[0] = DensName;
  DataLabel[1] = TEName;
  DataLabel[2] = Vel1Name;
  DataLabel[3] = Vel2Name;
  DataLabel[4] = Vel3Name;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ShockTubeDirection     = %"ISYM"\n", ShockTubeDirection);
    fprintf(Outfptr, "ShockTubeBoundary      = %"FSYM"\n\n", ShockTubeBoundary);
 
    fprintf(Outfptr, "ShockTubeLeftDensity   = %"FSYM"\n", ShockTubeDensity[0]);
    fprintf(Outfptr, "ShockTubeLeftPressure  = %"FSYM"\n", ShockTubePressure[0]);
    fprintf(Outfptr, "ShockTubeLeftVelocity  = %"FSYM"\n\n", ShockTubeVelocity[0]);
 
    fprintf(Outfptr, "ShockTubeRightDensity  = %"FSYM"\n", ShockTubeDensity[1]);
    fprintf(Outfptr, "ShockTubeRightPressure = %"FSYM"\n", ShockTubePressure[1]);
    fprintf(Outfptr, "ShockTubeRightVelocity = %"FSYM"\n\n", ShockTubeVelocity[1]);
  }
 
  return SUCCESS;
}
