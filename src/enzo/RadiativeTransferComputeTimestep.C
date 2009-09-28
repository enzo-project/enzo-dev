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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level)
{

  const int MaxStepsPerHydroStep = 8;
  const float PhotonCourantFactor = 1.0;

  LevelHierarchyEntry *Temp;
  bool InitialTimestep;
  int l, maxLevel;
  FLOAT OldTime, HydroTime;
  float ThisPhotonDT;
  bool FoundRadiation = false;

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

  dtPhoton = 10*huge_number;

  // Calculate timestep by limiting to a 50% max change in HII
  if (RadiativeTransferHIIRestrictedTimestep) {

    float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
      VelocityUnits = 1, TimeUnits = 1;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
	     &VelocityUnits, LevelArray[level]->GridData->ReturnTime());

    for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
      for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
	ThisPhotonDT = Temp->GridData->
	  ComputeRT_TimeStep2(DensityUnits, LengthUnits);
	dtPhoton = min(dtPhoton, ThisPhotonDT);
	if (Temp->GridData->RadiationPresent() == TRUE)
	  FoundRadiation = true;
      }
    CommunicationMinValue(dtPhoton);

  } // ENDIF

  // if not (FoundRadiation will always be false if the
  // RTUseHIIRestricted parameter isn't set) or no radiation, use
  // hydro timestep on finest level
  if (!FoundRadiation) {
    if (maxLevel > 0) {
      for (Temp = LevelArray[maxLevel]; Temp; Temp = Temp->NextGridThisLevel) {
	ThisPhotonDT = Temp->GridData->ComputeRT_TimeStep();
	dtPhoton = min(dtPhoton, ThisPhotonDT);
      } // ENDFOR grids
      dtPhoton = CommunicationMinValue(dtPhoton);
    } else
      dtPhoton = dtLevelAbove;
    dtPhoton *= PhotonCourantFactor;

    // Ensure that not too many photon timesteps are taken per hydro step
    HydroTime = LevelArray[maxLevel]->GridData->ReturnTime();
    dtPhoton = max(dtPhoton, 
		   (HydroTime - PhotonTime) / MaxStepsPerHydroStep);

  } // ENDELSE

  /* Do not go past the level-0 time (+timestep) */

  HydroTime = LevelArray[0]->GridData->ReturnTime();
  if (level == 0)
    HydroTime += LevelArray[0]->GridData->ReturnTimeStep();
  dtPhoton = min(HydroTime - PhotonTime, dtPhoton);

  if (InitialTimestep && !MetaData->FirstTimestepAfterRestart)
    dtPhoton = min(dtPhoton, dtLevelAbove);
  
  return SUCCESS;

}
