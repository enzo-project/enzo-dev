/***********************************************************************
/
/  ACTIVE PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- Converted to active particles.
/
/  PURPOSE: Prepare the active particles before going through the main
/  EvolveLevel grid loop.
/
************************************************************************/

#include "preincludes.h"
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
#include "ActiveParticle.h"

#define NO_DEBUG

int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[]);

#ifdef TRANSFER
void DeleteGlobalRadiationSources(void);
#endif

int ActiveParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			     int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			     int ThisLevel)
{

  int i;

  /* Return if this does not concern us */
  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleInitialize");

  int *TotalActiveParticleCount = new int[NumberOfGrids]();

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles;// - NumberOfActiveParticles;

  if (NextActiveParticleID == INT_UNDEFINED)
    NextActiveParticleID = NumberOfOtherParticles + NumberOfActiveParticles;

  bool CallEvolvePhotons = false;

#ifdef TRANSFER
  /* If radiation sources exist, delete them */
  if (RadiativeTransfer == TRUE) {

    DeleteGlobalRadiationSources();

    GlobalRadiationSources = new RadiationSourceEntry;
    GlobalRadiationSources->NextSource = NULL;
    GlobalRadiationSources->PreviousSource = NULL;

    /* Determine whether EvolvePhotons will be called this timestep */

    FLOAT GridTime = LevelArray[ThisLevel]->GridData->ReturnTime();
    float dt = LevelArray[ThisLevel]->GridData->ReturnTimeStep();
    if ((GridTime+dt > PhotonTime && LevelArray[ThisLevel+1] == NULL)
	|| MetaData->FirstTimestepAfterRestart)
      CallEvolvePhotons = true;
  }
#endif /* TRANSFER */

  /* Call initialization routines for each active particle type */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();

    ActiveParticleTypeToEvaluate->BeforeEvolveLevel(Grids, MetaData, NumberOfGrids, LevelArray, 
						    ThisLevel, CallEvolvePhotons, 
						    TotalActiveParticleCount,
						    ActiveParticleID);

  }

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

  delete [] TotalActiveParticleCount;

  LCAPERF_STOP("ActiveParticleInitialize");
  return SUCCESS;

}
