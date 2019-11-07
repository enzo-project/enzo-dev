/***********************************************************************
/
/  DISABLE THE ASSOCIATED ACTIVE PARTICLE (Delete frrom particle list)
/
/  written by: John Wise
/  date:       December, 2009
/  modified1:  Nathan Goldbaum, December 2011 (porting to active particles)
/
************************************************************************/

#include "preincludes.h"

#ifdef USE_MPI
#endif
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
#include "ActiveParticle.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int ActiveParticleType::DisableParticle(LevelHierarchyEntry *LevelArray[], int NewProcessorNumber)
{

  int i, ID, nPart, NumberOfGrids, changedGrid = INT_UNDEFINED, foundAP = FALSE,
    foundP = FALSE;
  HierarchyEntry **Grids;
  
  NumberOfGrids = GenerateGridArray(LevelArray, this->level, &Grids);
  for (i = 0; i < NumberOfGrids; i++) {
    ID = this->ReturnID();
    foundAP = Grids[i]->GridData->RemoveActiveParticle(ID, NewProcessorNumber);
    if (foundAP) {
      changedGrid = i;
      break;
    }
  } // ENDFOR grids

  delete [] Grids;

  return SUCCESS;

}
