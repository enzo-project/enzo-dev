/***********************************************************************
/
/  ACTIVE PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- converting to active particles
/
/  PURPOSE: Contains all routines to finalize the active particles.
/
************************************************************************/
#ifdef USE_MPI
#endif
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
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

#define NO_DEBUG

/* prototypes */

int CommunicationUpdateActiveParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids,
					 int TotalActiveParticleCountPrevious[]);

int ActiveParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int level, int NumberOfNewActiveParticles[])
{
  int i;

  if (EnabledActiveParticlesCount == 0) return SUCCESS;
  FLOAT TimeNow = LevelArray[level]->GridData->ReturnTime();
  float Timestep = LevelArray[level]->GridData->ReturnTimeStep();

  LCAPERF_START("ActiveParticleFinalize");

  /* Update the active particle counters. */

  CommunicationUpdateActiveParticleCount(Grids, MetaData, NumberOfGrids,
					 NumberOfNewActiveParticles);

#ifdef DEBUG
  
  int nParticles;
  ActiveParticleList ParticleList;

  ActiveParticleFindAll(LevelArray, &nParticles, 0, ParticleList);

  if (nParticles > 0) {
    PINT IDList[nParticles];
    for (i = 0; i < nParticles; i++)
      IDList[i] = ParticleList[i]->ReturnID();
    std::sort(IDList, IDList + sizeof(IDList)/sizeof(IDList[0]));
    for (i = 0; i < nParticles-1; i++)
      if (IDList[i] == IDList[i+1]) {
	ENZO_FAIL("Two active particles have identical IDs"); }
  }

#endif

  /* Call finalization routines for each active particle type  */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    

    ActiveParticleTypeToEvaluate->
      AfterEvolveLevel(Grids,MetaData,NumberOfGrids,LevelArray, 
		       level, NumberOfNewActiveParticles, ActiveParticleID);

  }

  LCAPERF_STOP("ActiveParticleFinalize");
  return SUCCESS;

}
