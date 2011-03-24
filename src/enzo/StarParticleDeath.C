/***********************************************************************
/
/  STAR PARTICLE DEATHS
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
/
/  PURPOSE: Kills any star particles and regular particles when a 
/           star dies.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include "performance.h"
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

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

Star *PopStar(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);
void DeleteStar(Star * &Node);

  /* Check for any stellar deaths */

int StarParticleDeath(LevelHierarchyEntry *LevelArray[], int level,
		      Star *&AllStars)
{

  int death;
  FLOAT TimeNow;
  Star *ThisStar, *MoveStar, *LastStar;

  LCAPERF_START("StarParticleDeath");
  ThisStar = AllStars;
  AllStars = NULL;
  LastStar = NULL;
  while (ThisStar) {
    TimeNow = LevelArray[ThisStar->ReturnLevel()]->GridData->ReturnTime();
    //TimeNow = LevelArray[level]->GridData->ReturnTime();
    death = ThisStar->HitEndpoint(TimeNow);
    MoveStar = PopStar(ThisStar);
    if (death == KILL_STAR) {
      MoveStar->MirrorToParticle();
      MoveStar->DeleteCopyInGridGlobal(LevelArray);
      DeleteStar(MoveStar);
    } else if (death == KILL_ALL) {
      // Never should be done.  Deleting particles messing the star
      // particle counts up.
      MoveStar->DeleteCopyInGrid();
      MoveStar->DeleteParticle(LevelArray);
      DeleteStar(MoveStar);
    } else {
      // Re-insert at the end of the list to keep the ordering the
      // same (for AddedFeedback array)
      if (LastStar == NULL)
	InsertStarAfter(AllStars, MoveStar);
      else
	InsertStarAfter(LastStar, MoveStar);
      LastStar = MoveStar;
    }

  } // ENDWHILE stars

  LCAPERF_STOP("StarParticleDeath");
  return SUCCESS;

}
