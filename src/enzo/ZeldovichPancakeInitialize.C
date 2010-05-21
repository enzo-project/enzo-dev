/***********************************************************************
/
/  INITIALIZE A ZELDOVICH PANCAKE
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
 
int ZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
 
  /* declarations */
 
  char line[MAX_LINE_LENGTH];
  int ret;
 
  /* Error check. */
 
  if (!ComovingCoordinates) {
    ENZO_FAIL("ComovingCoordinates must be TRUE!\n");
  }
 
  if (!SelfGravity)
    fprintf(stderr, "ZeldovichPancake: gravity is off!?!\n");
  if (CellFlaggingMethod[0] < 2)
    fprintf(stderr, "ZeldovichPancake: check CellFlaggingMethod.\n");
 
  /* set default parameters */
 
  int   ZeldovichPancakeDirection          = 0;    // along the x-axis
  float ZeldovichPancakeCentralOffset      = 0.0;  // no offset
  float ZeldovichPancakeOmegaBaryonNow     = 1.0;  // standard
  float ZeldovichPancakeOmegaCDMNow        = 0.0;  // no dark matter
  float ZeldovichPancakeCollapseRedshift   = 1.0;  // free parameter
  float ZeldovichPancakeInitialTemperature = 100;  // whatever
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "ZeldovichPancakeDirection = %"ISYM,
		  &ZeldovichPancakeDirection);
    ret += sscanf(line, "ZeldovichPancakeCentralOffset = %"FSYM,
		  &ZeldovichPancakeCentralOffset);
    ret += sscanf(line, "ZeldovichPancakeOmegaBaryonNow = %"FSYM,
		  &ZeldovichPancakeOmegaBaryonNow);
    ret += sscanf(line, "ZeldovichPancakeOmegaCDMNow = %"FSYM,
		  &ZeldovichPancakeOmegaCDMNow);
    ret += sscanf(line, "ZeldovichPancakeCollapseRedshift = %"FSYM,
		  &ZeldovichPancakeCollapseRedshift);
    ret += sscanf(line, "ZeldovichPancakeInitialTemperature = %"FSYM,
		  &ZeldovichPancakeInitialTemperature);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "ZeldovichPancake"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* set up grid */
 
  if (TopGrid.GridData->ZeldovichPancakeInitializeGrid(
					  ZeldovichPancakeDirection,
					  ZeldovichPancakeCentralOffset,
					  ZeldovichPancakeOmegaBaryonNow,
					  ZeldovichPancakeOmegaCDMNow,
					  ZeldovichPancakeCollapseRedshift,
					  ZeldovichPancakeInitialTemperature
						       ) == FAIL) {
    ENZO_FAIL("Error in ZeldovichPancakeInitializeGrid.\n");
  }
 
  /* set up field names and units */
 
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
 
  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "ZeldovichPancakeDirection          = %"ISYM"\n",
	    ZeldovichPancakeDirection);
    fprintf(Outfptr, "ZeldovichPancakeCentralOffset      = %"FSYM"\n",
	    ZeldovichPancakeCentralOffset);
    fprintf(Outfptr, "ZeldovichPancakeOmegaBaryonNow     = %"FSYM"\n",
	    ZeldovichPancakeOmegaBaryonNow);
    fprintf(Outfptr, "ZeldovichPancakeOmegaCDMNow        = %"FSYM"\n",
	    ZeldovichPancakeOmegaCDMNow);
    fprintf(Outfptr, "ZeldovichPancakeCollapseRedshift   = %"FSYM"\n",
	    ZeldovichPancakeCollapseRedshift);
    fprintf(Outfptr, "ZeldovichPancakeInitialTemperature = %"FSYM"\n\n",
	    ZeldovichPancakeInitialTemperature);
  }
 
  return SUCCESS;
}
