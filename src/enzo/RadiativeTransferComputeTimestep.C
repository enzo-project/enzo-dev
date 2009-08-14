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

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "CommunicationUtilities.h"

extern int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove)
{

  const int MaxStepsPerHydroStep = 8;

  LevelHierarchyEntry *Temp;
  bool InitialTimestep;
  int l, level, maxLevel;
  FLOAT HydroTime;
  float ThisPhotonDT;

  // Search for the maximum level
  for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
    if (LevelArray[l] == NULL) {
      maxLevel = l-1;
      break;
    }

  // Determine if this is the first timestep (not in restart)
  InitialTimestep = true;
  for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
    if (LevelCycleCount[l] > 0) {
      InitialTimestep = false;
      break;
    }

  // Now calculate the timestep on the max. level
  dtPhoton = 1e20;
  for (Temp = LevelArray[maxLevel]; Temp; Temp = Temp->NextGridThisLevel) {
    ThisPhotonDT = Temp->GridData->ComputeTimeStep();
    dtPhoton = min(dtPhoton, ThisPhotonDT);
  } // ENDFOR grids
  dtPhoton = CommunicationMinValue(dtPhoton);

  // Ensure that not too many photon timesteps are taken per hydro step
  HydroTime = LevelArray[maxLevel]->GridData->ReturnTime();
  dtPhoton = max(dtPhoton, 
		 (HydroTime - PhotonTime) / MaxStepsPerHydroStep);
  if (InitialTimestep && !MetaData->FirstTimestepAfterRestart)
    dtPhoton = min(dtPhoton, dtLevelAbove);
  
  return SUCCESS;

}
