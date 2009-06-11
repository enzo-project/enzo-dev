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
#include "LevelHierarchy.h"
 
int ProcMapper(LevelHierarchyEntry *LevelArray[])
{
 
  int proc;
  int level;
  int gridcounter;
 
  gridcounter = 0;
 
  // Walk the grids
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* Loop over all the grids. */
 
    LevelHierarchyEntry *Temp = LevelArray[level];
 
    while (Temp != NULL) {
 
      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
 
      while (Temp2 != NULL) {
        proc = Temp2->GridData->ReturnProcessorNumber();
        fprintf(stderr, "Calling proc %"ISYM" : this proc %"ISYM"\n", MyProcessorNumber, proc);
 
        Temp2 = Temp2->NextGridThisLevel;
      }
 
      /* Next grid on this level. */
 
      Temp = Temp->NextGridThisLevel;
 
    } // end loop over grids
 
  } // end loop over levels
 
  return SUCCESS;
 
}
