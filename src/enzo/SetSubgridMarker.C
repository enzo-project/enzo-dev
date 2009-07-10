/***********************************************************************
/
/  SetSubgridMarker FUNCTION
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:  
/
/  PURPOSE: Set the SubgridMarker field for all grids on finer 
/           levels than this one
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level)
{
  int i, grid;
  LevelHierarchyEntry *Temp;
  HierarchyEntry *Subgrid;
  int NumberOfGrids;

  for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++)  {

    Temp = LevelArray[i];
    while (Temp != NULL) {

      // first the grid marks itself
      Temp->GridData->
	SetSubgridMarkerFromSubgrid(Temp->GridData, Temp->GridData);
      
      Subgrid = Temp->GridHierarchyEntry->NextGridNextLevel;

      while (Subgrid != NULL) {
	Temp->GridData->
	  SetSubgridMarkerFromSubgrid(Subgrid->GridData, Temp->GridData);
	Subgrid = Subgrid->NextGridThisLevel;
      } // ENDWHILE Subgrid

      Temp = Temp->NextGridThisLevel;

    } // ENDWHILE grids
  } // ENDFOR levels

  return SUCCESS;
}
