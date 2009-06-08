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

#ifndef LEVEL_HIERARCHY_ENTRY_DEFINED__
#define LEVEL_HIERARCHY_ENTRY_DEFINED__

#include "Grid.h"

struct LevelHierarchyEntry
{
  LevelHierarchyEntry *NextGridThisLevel;  /* next entry on this level */
  grid                *GridData;           /* pointer to this entry's grid */
  HierarchyEntry      *GridHierarchyEntry; /* pointer into hierarchy */
};

#endif
