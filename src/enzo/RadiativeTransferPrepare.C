/***********************************************************************
/
/  PREPARE THE RADIATIVE TRANSFER MODULE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/ PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "Star.h"

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove);
int StarParticleRadTransfer(LevelHierarchyEntry *LevelArray[], int level,
			    Star *AllStars);
int RestartPhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		   Star *AllStars);

int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove)
{

  /* Determine the photon timestep */

  if (RadiativeTransferComputeTimestep(LevelArray, MetaData, 
				       dtLevelAbove) == FAIL) {
    fprintf(stderr, "Error in RadiativeTransferComputeTimestep.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  /* Convert star particles into radiation sources only if we're going
     into EvolvePhotons */

  if (dtPhoton > 0.0 && LevelArray[level]->GridData->ReturnTime() >= PhotonTime)
    if (StarParticleRadTransfer(LevelArray, level, AllStars) == FAIL) {
      fprintf(stderr, "Error in StarParticleRadTransfer.\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
  
  /* If the first timestep after restart and we have radiation
     sources, go back a light crossing time of the box and run
     EvolvePhotons to populate the grids with the proper rates. */

  if (MetaData->FirstTimestepAfterRestart == TRUE)
    if (RestartPhotons(MetaData, LevelArray, AllStars) == FAIL) {
      fprintf(stderr, "Error in RestartPhotons.\n");
      ENZO_FAIL("Error in: "__FILE__);
    }


  return SUCCESS;

}
