/***********************************************************************
/
/  CONVERT DATA HIERARCHY TO NUMPY FORMAT
/
/  written by: Matthew Turk
/  date:       September, 2008
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function writes out the data hierarchy (TopGrid)

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


#ifdef USE_PYTHON

int ExposeDataHierarchy(TopGridData *MetaData, HierarchyEntry *Grid, 
		       int &GridID, FLOAT WriteTime, int reset, int ParentID, int level)
{

  int OriginalID, NextGridThisLevelID, NextGridNextLevelID;
  int flagged, noParent = 0;
  static PyArrayObject *container[11];
  
  /* Borrowed reference */
  if(reset == 1) {
        int i = 0;
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridDimensions");
        container[i++] =  (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridStartIndices");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridEndIndices");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridLeftEdge");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridRightEdge");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridLevels");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridTimes");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridOldTimes");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridProcs");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridNumberOfParticles");
        container[i++] = (PyArrayObject *)
            PyDict_GetItemString(hierarchy_information, "GridParentIDs");
      }

  /* Write out header info for this grid */

  OriginalID = GridID;

  /* We should have a call for interpolation here */

  Grid->GridData->ConvertToNumpy(GridID, container, ParentID, level, WriteTime);

  /* Write out pointer information for the next grid this level */

  NextGridThisLevelID = GridID + 1;
  if (Grid->NextGridThisLevel == NULL) NextGridThisLevelID = 0;

  if (NextGridThisLevelID != 0) {
    GridID++;
    if (ExposeDataHierarchy(MetaData, Grid->NextGridThisLevel, GridID,
			   WriteTime, 0, ParentID, level) == FAIL) {
      fprintf(stderr, "Error in ExposeDataHierarchy(1).\n");
      return FAIL;
    }
  }

  /* Write out pointer information for the next grid next level */

  NextGridNextLevelID = GridID + 1;
  if (Grid->NextGridNextLevel == NULL) NextGridNextLevelID = 0;

  if (NextGridNextLevelID != 0) {
    GridID++;
    if (ExposeDataHierarchy(MetaData, Grid->NextGridNextLevel, GridID,
			   WriteTime, 0, OriginalID, level+1) == FAIL) {
      fprintf(stderr, "Error in ExposeDataHierarchy(1).\n");
      return FAIL;
    }
  }

  return SUCCESS;
}
#endif
