/***********************************************************************
/
/  INITIALIZE A 1D PRESSURELESS COLLAPSE
/
/  written by: Greg Bryan
/  date:       April, 1995
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
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
 
int PressurelessCollapseInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int ret;
 
  /* Error check. */
 
  if (ComovingCoordinates)
    fprintf(stderr, "PressurelessCollapse: ComovingCoordinates are ON?\n");
  if (!SelfGravity)
    fprintf(stderr, "PressurelessCollapse: gravity is off!?!");
  if (CellFlaggingMethod[0] < 2)
    fprintf(stderr, "PressurelessCollapse: check CellFlaggingMethod.\n");
  if (MetaData.GravityBoundary != TopGridIsolated)
    fprintf(stderr, "PressurelessCollapse: check GravityBoundary.\n");
  if (PressureFree != TRUE)
    fprintf(stderr, "PressurelessCollapse: PressureFree is not ON!\n");
 
  /* set default parameters */
 
  int   PressurelessCollapseDirection          = 0;    // along the x-axis
  float PressurelessCollapseInitialDensity     = 1;
  int   PressurelessCollapseNumberOfCells      = INT_UNDEFINED;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "PressurelessCollapseDirection = %"ISYM,
		  &PressurelessCollapseDirection);
    ret += sscanf(line, "PressurelessCollapseInitialDensity = %"FSYM,
		  &PressurelessCollapseInitialDensity);
    ret += sscanf(line, "PressurelessCollapseNumberOfCells = %"ISYM,
		  &PressurelessCollapseNumberOfCells);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "PressurelessCollapse"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */
 
  if (TopGrid.GridData->PressurelessCollapseInitializeGrid(
					  PressurelessCollapseDirection,
					  PressurelessCollapseInitialDensity,
					  PressurelessCollapseNumberOfCells
						       ) == FAIL) {
    ENZO_FAIL("Error in PressurelessCollapseInitializeGrid.\n");
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

    fprintf(Outfptr, "PressurelessCollapseDirection      = %"ISYM"\n",
	    PressurelessCollapseDirection);
    fprintf(Outfptr, "PressurelessCollapseInitialDensity = %"FSYM"\n\n",
	    PressurelessCollapseInitialDensity);
  }
 
  return SUCCESS;
}
