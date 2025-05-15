/***********************************************************************
/
/  CONVERT DATA HIERARCHY TO LIBYT FORMAT
/
/  written by: Matthew Turk
/  date:       April, 2023
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function writes out the data hierarchy (TopGrid)
#ifdef USE_LIBYT
#include "libyt.h"
#endif

#include <stdlib.h>
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

/* function prototypes */


#ifdef USE_LIBYT

int ExposeHierarchyToLibyt(TopGridData *MetaData, HierarchyEntry *Grid, 
               int &GridID, int &LocalGridID, FLOAT WriteTime, int ParentID,
               int level, yt_grid *GridInfoArray)
{

  int OriginalID;

  /* Write out header info for this grid */

  OriginalID = GridID;

  /* We should have a call for interpolation here */

  if (Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      /* We do a similar */
    Grid->GridData->ConvertToLibyt(LocalGridID, GridID, ParentID, level,
            GridInfoArray[LocalGridID-1]);
    LocalGridID++;
  }

  if (Grid->NextGridThisLevel != NULL) {
    GridID++;
    if (ExposeHierarchyToLibyt(MetaData, Grid->NextGridThisLevel, GridID, LocalGridID,
			   WriteTime, ParentID, level, GridInfoArray) == FAIL) {
      fprintf(stderr, "Error in ExposeHierarchyToLibyt(1).\n");
      return FAIL;
    }
  }

  if (Grid->NextGridNextLevel != NULL) {
    GridID++;
    if (ExposeHierarchyToLibyt(MetaData, Grid->NextGridNextLevel, GridID, LocalGridID,
			   WriteTime, OriginalID, level+1, GridInfoArray) == FAIL) {
      fprintf(stderr, "Error in ExposeDataHierarchy(2).\n");
      return FAIL;
    }
  }

  return SUCCESS;
}
#endif
