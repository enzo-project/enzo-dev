/***********************************************************************
/
/  REMOVE THE STAR PARTICLE AND ADJUST NUMBEROFSTARS IN ALL PROCESSORS
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
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
#include "CommunicationUtilities.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int Star::DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[])
{

  int i, nPart, n, NumberOfGrids, changedGrid = INT_UNDEFINED, found = FALSE;
  HierarchyEntry **Grids;
  
  NumberOfGrids = GenerateGridArray(LevelArray, this->level, &Grids);
  if (this->CurrentGrid != NULL)  // i.e. star on local processor
    for (i = 0; i < NumberOfGrids; i++)
      if (Grids[i]->GridData == this->CurrentGrid) {
	changedGrid = i;
	break;
      }

#ifdef USE_MPI
  CommunicationAllReduceValues(&changedGrid, 1, MPI_MAX);
#endif

  if (changedGrid == INT_UNDEFINED) {
    if (debug)
      fprintf(stdout, "Star::DeleteCopyInGridGlobal: WARNING -- "
	      "star %"ISYM" not found...\n", this->Identifier);
    delete [] Grids;
    return SUCCESS;
  }

  if (this->CurrentGrid != NULL) // meaning that the grid is on this processor
    this->DeleteCopyInGrid();
  Grids[changedGrid]->GridData->NumberOfStars--;

  delete [] Grids;

  return SUCCESS;

}
