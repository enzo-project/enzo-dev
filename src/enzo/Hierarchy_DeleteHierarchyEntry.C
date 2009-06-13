/***********************************************************************
/
/  HIERARCHY STRUCTURE (DELETE HIERARCHY ENTRY & ALL ITS CHILDREN & SIBLINGS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// This function recursively deletes a point in the hierarchy and all
//   it's children and siblings (but not the grid data itself, so be careful).

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
 
void DeleteHierarchyEntry(HierarchyEntry *Grid)
{
 
  /* Delete siblings (i.e. grids on the same level). */
 
  if (Grid->NextGridThisLevel != NULL)
    DeleteHierarchyEntry(Grid->NextGridThisLevel);
 
  /* Delete children. */
 
  if (Grid->NextGridNextLevel != NULL)
    DeleteHierarchyEntry(Grid->NextGridNextLevel);
 
  /* Delete this HierarchyEntry (but not the grid data itself). */
 
  delete Grid;
}
