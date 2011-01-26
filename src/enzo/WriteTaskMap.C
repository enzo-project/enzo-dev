/***********************************************************************
/
/  OUTPUT GRID-TO-TASK MEMORY MAP
/
/  written by: Robert Harkness
/  date:       January 2006
/  modified1:
/
/  PURPOSE: DO NOT USE if GRIDS:Processors != 1:1
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
 
 
 
 
int WriteTaskMap(FILE *fptr, HierarchyEntry *Grid,
		 char* base_name, int &GridID, FLOAT WriteTime)
{
 
  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
 
  OriginalID = GridID;
 
  if (Grid->GridData->WriteTaskMap(fptr, base_name, GridID) == FAIL) {
    ENZO_FAIL("Error in grid->WriteTaskMap.\n");
  }

  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;

  if (NextGridThisLevelID != 0) {
    GridID++;
    if (WriteTaskMap(fptr, Grid->NextGridThisLevel, base_name, GridID, WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteTaskMap(1).\n");
    }
  }
 
  NextGridNextLevelID = GridID + 1;
  if (Grid->NextGridNextLevel == NULL) NextGridNextLevelID = 0;
 
  if (NextGridNextLevelID != 0) {
    GridID++;
    if (WriteTaskMap(fptr, Grid->NextGridNextLevel, base_name, GridID, WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteTaskMap (1).\n");

    }
  }
 
  return SUCCESS;
}
