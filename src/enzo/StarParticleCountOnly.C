/***********************************************************************
/
/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: First synchronizes particle information in the normal and 
/           star particles.  Then we make a global particle list, which
/           simplifies adding the feedback to the all of the grids.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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

int StarParticleCountOnly(LevelHierarchyEntry *LevelArray[])
{

  int level, LocalNumberOfStars;
  float LocalMinLifetime;
  LevelHierarchyEntry *Temp;

  if (StarParticleFeedback == 0 && StarParticleCreation == 0)
    return SUCCESS;

  LocalMinLifetime = 1e20;
  LocalNumberOfStars = 0;

  /* Find statistics on local grids */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->ReturnStarStatistics(LocalNumberOfStars, LocalMinLifetime);

  /* Reduce the statistics */

  CommunicationAllSumValues(&LocalNumberOfStars, 1);
#ifdef USE_MPI
  CommunicationAllReduceValues(&LocalMinLifetime, 1, MPI_MIN);
#endif

  // Store them in the global variables
  G_TotalNumberOfStars = LocalNumberOfStars;
  minStarLifetime = LocalMinLifetime;

  return SUCCESS;

}
