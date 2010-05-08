/***********************************************************************
/
/  REMOVE THE ASSOCIATED NORMAL PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: shouldn't be used any more in general, instead try 
/             Star_DisableParticle
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

int Star::DeleteParticle(LevelHierarchyEntry *LevelArray[])
{

  int i, nPart, NumberOfGrids, changedGrid = INT_UNDEFINED, found = FALSE;
  HierarchyEntry **Grids;
  
  NumberOfGrids = GenerateGridArray(LevelArray, this->level, &Grids);
  for (i = 0; i < NumberOfGrids; i++) {
    found = Grids[i]->GridData->RemoveParticle(this->Identifier);
    if (found) {
      changedGrid = i;
      break;
    }
  } // ENDFOR grids

  /* Now clean up deleted particle on the local processor and adjust
     NumberOfParticles on others */

#ifdef USE_MPI
  CommunicationAllReduceValues(&changedGrid, 1, MPI_MAX);
#endif

  if (changedGrid == INT_UNDEFINED) {
    if (debug)
      fprintf(stdout, "RemoveParticles: WARNING -- "
	      "particle %"ISYM" not found...\n", this->Identifier);
    delete [] Grids;
    return SUCCESS;
  }

  if (found) // meaning that the grid is on this processor
    Grids[changedGrid]->GridData->CleanUpMovedParticles();
  else
    Grids[changedGrid]->GridData->NumberOfParticles--;

  Grids[changedGrid]->GridData->NumberOfStars--;

  delete [] Grids;

  return SUCCESS;

}
