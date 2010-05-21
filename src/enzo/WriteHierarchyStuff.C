/***********************************************************************
/
/  WRITE OUT THE HIERARCHY
/
/  written by: Robert Harkness
/  date:       July, 2003
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// This function writes out the hierarchy (TopGrid)
 
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
 
 
 
 
int WriteHierarchyStuff(FILE *fptr, HierarchyEntry *Grid,
		        char* base_name, int &GridID, FLOAT WriteTime)
{
 
  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
 
  /* Write out header info for this grid */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "\nGrid = %"ISYM"\n", GridID);
 
  OriginalID = GridID;
 
  /* Write out grid data for this grid (if WriteTime is < 0 then output
     at the grid's time, otherwise interpolate to WriteTime and output). */
 
    if (Grid->GridData->WriteStuff(fptr, base_name, GridID) == FAIL) {
      ENZO_FAIL("Error in grid->WriteGrid.\n");
    }
 
  /* Write out pointer information for the next grid this level */
 
  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "Pointer: Grid[%"ISYM"]->NextGridThisLevel = %"ISYM"\n", OriginalID,
	    NextGridThisLevelID);
 
  if (NextGridThisLevelID != 0) {
    GridID++;
    if (WriteHierarchyStuff(fptr, Grid->NextGridThisLevel, base_name, GridID,
			   WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteHierarchyStuff(1).\n");
    }
  }
 
  /* Write out pointer information for the next grid next level */
 
  NextGridNextLevelID = GridID + 1;
  if (Grid->NextGridNextLevel == NULL) NextGridNextLevelID = 0;
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "Pointer: Grid[%"ISYM"]->NextGridNextLevel = %"ISYM"\n", OriginalID,
	    NextGridNextLevelID);
 
  if (NextGridNextLevelID != 0) {
    GridID++;
    if (WriteHierarchyStuff(fptr, Grid->NextGridNextLevel, base_name, GridID,
			   WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteHierarchyStuff(1).\n");

    }
  }
 
  return SUCCESS;
}
