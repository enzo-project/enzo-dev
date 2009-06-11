/***********************************************************************
/
/  STAR PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Contains all routines to initialize the star particles.
/
************************************************************************/

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

int StarParticleFindAll(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int StarParticleMergeNew(LevelHierarchyEntry *LevelArray[], Star *&AllStars);

int StarParticleInitialize(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			   TopGridData *MetaData, Star *&AllStars)
{

  /* Return if this does not concern us */
  if (!(StarParticleCreation || StarParticleFeedback)) return SUCCESS;

  int level;
  Star *cstar;
  LevelHierarchyEntry *Temp;
  FLOAT TimeNow = LevelArray[ThisLevel]->GridData->ReturnTime();

  /* Initialize all star particles if this is a restart */

  if (MetaData->FirstTimestepAfterRestart)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->FindAllStarParticles(level) == FAIL) {
	  fprintf(stderr, "Error in grid::FindAllStarParticles.\n");
	  ENZO_FAIL("");
	}

  /* Create a master list of all star particles */

  if (StarParticleFindAll(LevelArray, AllStars) == FAIL) {
    fprintf(stderr, "Error in StarParticleFindAll.\n");
    ENZO_FAIL("");
  }

  /* Merge any newly created, clustered particles */

  if (StarParticleMergeNew(LevelArray, AllStars) == FAIL) {
    fprintf(stderr, "Error in StarParticleMergeNew.\n");
    ENZO_FAIL("");
  }

  /* 
     Set feedback flags.  

     Sync all star and normal particles that are stored in the grids
     to the global list (AllStars) so these changes are reflected
     there.
  */

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {
    cstar->SetFeedbackFlag(TimeNow);
    cstar->CopyToGrid();
    cstar->MirrorToParticle();
  }

  return SUCCESS;

}
