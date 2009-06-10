#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "LevelHierarchy.h"

int ReassignSuperSources(LevelHierarchyEntry *LevelArray[])
{

  if (LevelArray == NULL)
    return SUCCESS;

  if (OldSourceClusteringTree == NULL)
    return SUCCESS;

  int level;
  LevelHierarchyEntry *Temp;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->ReassignSuperSources() == FAIL) {
	fprintf(stderr, "Error in grid::ReassignSuperSources.\n");
	ENZO_FAIL("Error in: "__FILE__);
      }

  return SUCCESS;

}
