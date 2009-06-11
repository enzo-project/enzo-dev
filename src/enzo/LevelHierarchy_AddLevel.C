/***********************************************************************
/
/  LEVEL HIERARCHY STRUCTURE AND ROUTINES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// This function adds the hierarchy entry Grid to the head of the linked
//   list begining at LevelArray[level].  It then recursively adds all
//   the subgrids of Grid.

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
#include "LevelHierarchy.h"
 
void AddLevel(LevelHierarchyEntry *LevelArray[], HierarchyEntry *Grid,
	      int level)
{
  LevelHierarchyEntry *ThisLevel;
 
  /* create a new LevelHierarchyEntry for the HierarchyEntry Grid
     and insert it into the head of the linked list (LevelArray[level]). */
 
  ThisLevel = new LevelHierarchyEntry;
  ThisLevel->GridData = Grid->GridData;
  ThisLevel->NextGridThisLevel = LevelArray[level];
  ThisLevel->GridHierarchyEntry = Grid;
  LevelArray[level] = ThisLevel;
 
  /* recursively call this for the next grid on this level. */
 
  if (Grid->NextGridThisLevel != NULL)
    AddLevel(LevelArray, Grid->NextGridThisLevel, level);
 
  /* ... and then descend the tree. */
 
  if (Grid->NextGridNextLevel != NULL)
    AddLevel(LevelArray, Grid->NextGridNextLevel, level+1);
}
