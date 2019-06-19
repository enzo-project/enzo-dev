/***********************************************************************
/
/  ACTIVE PARTICLE RESET ACCELERATION
/
/  written by: John Regan
/  date:       June, 2016
/  modified1:  
/
/  PURPOSE: This allows active particles to remain fixed in comoving space
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

#define DEBUG 0


/*
 * Reset the acceleration for the given active particle
 * What is passed in is a 1 D arrays referring to a single direction 
 * for all active particles. We then loop over all active particles 
 * and call the API hook. 
 */
int ActiveParticleResetAccelerations(float *ActiveParticleAcceleration)
{

  /* Call reset accelerations routines for each active particle type  */

  for (int i = 0 ; i < EnabledActiveParticlesCount; i++) {
    int ActiveParticleID;
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleID = ActiveParticleTypeToEvaluate->GetEnabledParticleID();
    
    /* 
     * For each enabled active particle reset the acceleration
     */
    ActiveParticleTypeToEvaluate->ResetAcceleration(&(ActiveParticleAcceleration[i]));
  }

  return SUCCESS;

}
