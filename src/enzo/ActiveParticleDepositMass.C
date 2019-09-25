/***********************************************************************
/
/  ACTIVE PARTICLE DEPOSIT MASS
/
/  written by: Stephen Skory
/  date:       October, 2012
/  modified1:  
/
/  PURPOSE: This allows active particles to deposit mass into the temporary
/           buffers for computing the gravitational potential.
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

#define NODEBUG

/* prototypes */

ActiveParticleList<ActiveParticleType> &ActiveParticleFindAll(
    LevelHierarchyEntry *LevelArray[], int *GlobalNumberOfActiveParticles, 
    int ActiveParticleIDToFind);

int ActiveParticleDepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int level)
{
  int i;

  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  LCAPERF_START("ActiveParticleFinalize");

#ifdef DEBUG
  
  int nParticles;
  ActiveParticleList<ActiveParticleType> ParticleList = 
    ActiveParticleFindAll(LevelArray, &nParticles, 0);

  if (nParticles > 0) {
    PINT IDList[nParticles];
    for (i = 0; i < nParticles; i++)
      IDList[i] = ParticleList[i]->ReturnID();
    std::sort(IDList, IDList + sizeof(IDList)/sizeof(IDList[0]));
    for (i = 0; i < nParticles-1; i++)
      if (IDList[i] == IDList[i+1]) {
	ENZO_FAIL("Two active particles have identical IDs"); }
  }

  if (NumberOfProcessors > 1)
    for (i = 0; i < nParticles; i++)
      delete ParticleList[i];

#endif

  /* Call mass deposit routines for each active particle type  */

  int ActiveParticleID;

  for (i = 0 ; i < EnabledActiveParticlesCount; i++) {
    
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    

    ActiveParticleTypeToEvaluate->
      DepositMass(Grids,MetaData,NumberOfGrids,LevelArray, 
		       level, ActiveParticleID);

  }

  LCAPERF_STOP("ActiveParticleFinalize");
  return SUCCESS;

}
