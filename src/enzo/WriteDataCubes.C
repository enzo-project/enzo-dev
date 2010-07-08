/***********************************************************************
/
/  UNIGRID CUBE ASSEMBLY
/
/  written by: Robert Harkness
/  date:       April 2004
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
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
 
 
 
 
int WriteDataCubes(HierarchyEntry *Grid, int TGdims[],
		   char* base_name, int &GridID, FLOAT WriteTime)
{
 
  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
 
  OriginalID = GridID;
 
  /* Write out grid data for this grid (if WriteTime is < 0 then output
     at the grid's time, otherwise interpolate to WriteTime and output). */
 
  if (WriteTime < 0) {
    if (Grid->GridData->WriteCube(base_name, GridID, TGdims) == FAIL) {
      ENZO_FAIL("Error in grid->WriteCube.\n");
    }
  }
  else
    if (Grid->GridData->WriteCubeInterpolate(WriteTime, base_name, GridID, TGdims) == FAIL) {
      ENZO_FAIL("Error in grid->WriteCubeInterpolate.\n");
    }
 
  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;
 
  if (NextGridThisLevelID != 0) {
    GridID++;
    if (WriteDataCubes(Grid->NextGridThisLevel, TGdims, base_name, GridID, WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteDataCubes (1).\n");

    }
  }
 
  return SUCCESS;
}
