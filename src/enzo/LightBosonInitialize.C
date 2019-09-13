/***********************************************************************
/
/  INITIALIZE AN FDM SIMULATION FOR A 1D TEST
/
/  written by: Xinyu Li
/  date:       
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
  char *RePsiName = "RePsi";
  char *ImPsiName = "ImPsi";
  char *FDMDensityName = "FDMDensity";

  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int dim, ret;
 
  /* set default parameters */
 
  FLOAT LightBosonCenter = 0.5;     //
  int LightBosonProblemType = 1;    
  
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "LightBosonCenter = %"PSYM, &LightBosonCenter);
    ret += sscanf(line, "LightBosonProblemType = %"ISYM, &LightBosonProblemType);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "LightBoson"))
      fprintf(stderr, "LightBoson warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */

  if (TopGrid.GridData->
      LightBosonInitializeGrid(LightBosonCenter, LightBosonProblemType) == FAIL) {
    ENZO_FAIL("Error in LightBosonInitializeGrid (called from LightBosonInitialize).\n");
  } 

    /* set up field names and units */

  DataLabel[0] = DensName;
  DataLabel[1] = Vel1Name;
  DataLabel[2] = Vel2Name;
  DataLabel[3] = Vel3Name;
  DataLabel[4] = TEName;
  DataLabel[5] = RePsiName;
  DataLabel[6] = ImPsiName;
  DataLabel[7] = FDMDensityName;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  DataUnits[5] = NULL;
  DataUnits[6] = NULL;
  DataUnits[7] = NULL;

  
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "LightBosonCenter     = %"PSYM"\n", LightBosonCenter);
    fprintf(Outfptr, "LightBosonProblemType = %"ISYM"\n", LightBosonProblemType);
  }

  return SUCCESS;
}
