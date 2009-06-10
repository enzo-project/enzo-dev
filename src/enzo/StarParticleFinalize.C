/***********************************************************************
/
/  STAR PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Contains all routines to finalize the star particles.
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

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int StarParticleAddFeedback(TopGridData *MetaData, 
			    LevelHierarchyEntry *LevelArray[], int level, 
			    Star *&AllStars);
int StarParticleAccretion(Star *&AllStars);
int StarParticleDeath(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
void DeleteStarList(Star * &Node);
#ifdef USE_MPI
int CommunicationReduceValues(float *Values, int Number, MPI_Op ReduceOperation);
#endif

int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars)
{

  if (!StarParticleCreation && !StarParticleFeedback)
    return SUCCESS;

  int l;
  float TotalMass;
  Star *ThisStar, *MoveStar;
  LevelHierarchyEntry *Temp;
  FLOAT TimeNow;

  /* Update the star particle counters. */

  if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					   NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in CommunicationUpdateStarParticleCount.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  /* Update position and velocity of star particles from the actual
     particles */

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar)
    ThisStar->UpdatePositionVelocity();

  /* Apply any stellar feedback onto the grids and add any gas to the
     accretion rates of the star particles */

  if (StarParticleAddFeedback(MetaData, LevelArray, level, 
			      AllStars) == FAIL) {
    fprintf(stderr, "Error in StarParticleAddFeedback.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  /* Update star particles for any accretion */

  if (StarParticleAccretion(AllStars) == FAIL) {
    fprintf(stderr, "Error in StarParticleAccretion.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  /* Collect all sink particles and report the total mass to STDOUT */
  
  if (STARMAKE_METHOD(SINK_PARTICLE) && level == MaximumRefinementLevel) {
    TotalMass = 0.0;
    for (l = 0; l <= MaximumRefinementLevel; l++)
      for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel)
	TotalMass += Temp->GridData->ReturnTotalSinkMass();
#ifdef USE_MPI
    CommunicationReduceValues(&TotalMass, 1, MPI_SUM);
#endif
    if (debug)
      fprintf(stdout, "SinkParticle: Time = %"GOUTSYM", TotalMass = %"GSYM"\n", 
	      TimeNow, TotalMass);
  }

  /* Check for any stellar deaths */

  if (StarParticleDeath(LevelArray, AllStars) == FAIL) {
    fprintf(stderr, "Error in StarParticleDeath.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  /* 
     If the new particles are above a specified mass threshold,
     "activate" them.  Then check for any stellar deaths.

     Sync all star and normal particles that are stored in the grids
     to the global list (AllStars) so these changes are reflected
     there. 
  */

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {
    TimeNow = LevelArray[ThisStar->ReturnLevel()]->GridData->ReturnTime();
    ThisStar->ActivateNewStar(TimeNow);
    ThisStar->ResetAccretion();
    ThisStar->CopyToGrid();
    ThisStar->MirrorToParticle();
  } // ENDFOR stars


  /* Delete the global star particle list, AllStars */

  DeleteStarList(AllStars);

  return SUCCESS;

}
