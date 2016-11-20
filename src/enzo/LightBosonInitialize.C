/***********************************************************************
/
/  INITIALIZE A SHOCK IN A BOX
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
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
 
int LightBosonInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";




  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, ret;
 
  /* set default parameters */
 
  FLOAT LightBosonCenter = 0.5;     //
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "LightBosonCenter = %"PSYM, &LightBosonCenter);
 
    /*ret += sscanf(line, "ShockInABoxLeftDensity = %"ESYM,
		  &ShockInABoxDensity[0]);
    ret += sscanf(line, "ShockInABoxLeftPressure = %"ESYM,
		  &ShockInABoxPressure[0]);
    ret += sscanf(line, "ShockInABoxLeftVelocity = %"ESYM,
		  &ShockInABoxVelocity[0]);
 
    ret += sscanf(line, "ShockInABoxRightDensity = %"ESYM,
		  &ShockInABoxDensity[1]);
    ret += sscanf(line, "ShockInABoxRightPressure = %"ESYM,
		  &ShockInABoxPressure[1]);
    ret += sscanf(line, "ShockInABoxRightVelocity = %"ESYM,
		  &ShockInABoxVelocity[1]);
 
    ret += sscanf(line, "ShockInABoxSubgridLeft = %"PSYM,
		  &ShockInABoxSubgridLeft);
    ret += sscanf(line, "ShockInABoxSubgridRight = %"PSYM,
		  &ShockInABoxSubgridRight);*/
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "LightBoson"))
      fprintf(stderr, "LightBoson warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */

  if (TopGrid.GridData->
      LightBosonInitializeGrid(LightBosonCenter) == FAIL) {
    ENZO_FAIL("Error in LightBosonInitializeGrid (called from LightBosonInitialize).\n");
  } 

    /* set up field names and units */


  DataLabel[0] = DensName;
  DataLabel[1] = Vel1Name;
  DataLabel[2] = Vel2Name;
  DataLabel[3] = Vel3Name;
  DataLabel[4] = TEName;


 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;


  
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "LightBosonCenter     = %"PSYM"\n", LightBosonCenter);
    /*fprintf(Outfptr, "ShockInABoxBoundary      = %"GOUTSYM"\n\n",
	    ShockInABoxBoundary);
 
    fprintf(Outfptr, "ShockInABoxLeftDensity   = %"ESYM"\n", ShockInABoxDensity[0]);
    fprintf(Outfptr, "ShockInABoxLeftPressure  = %"ESYM"\n",
	    ShockInABoxPressure[0]);
    fprintf(Outfptr, "ShockInABoxLeftVelocity  = %"ESYM"\n\n",
	    ShockInABoxVelocity[0]);
 
    fprintf(Outfptr, "ShockInABoxRightDensity  = %"ESYM"\n", ShockInABoxDensity[1]);
    fprintf(Outfptr, "ShockInABoxRightPressure = %"ESYM"\n",
	    ShockInABoxPressure[1]);
    fprintf(Outfptr, "ShockInABoxRightVelocity = %"ESYM"\n\n",
	    ShockInABoxVelocity[1]);*/
  }
 
  return SUCCESS;
}
