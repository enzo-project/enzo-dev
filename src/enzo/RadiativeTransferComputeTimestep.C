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
//int LastTimestepUseHII = FALSE;
float LastPhotonDT = -1;

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level)
{

  if (!RadiativeTransferAdaptiveTimestep && dtPhoton != FLOAT_UNDEFINED)
    return SUCCESS;

  const int MaxStepsPerHydroStep = 8;
  const float PhotonCourantFactor = 1.0;

  // Restrict the increase in dtPhoton to this factor
  const float MaxDTChange = 1.0;

  LevelHierarchyEntry *Temp;
  bool InitialTimestep;
  int l, maxLevel;
  FLOAT OldTime, HydroTime;
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

//  LastPhotonDT = (MetaData->FirstTimestepAfterRestart && 
//		  LastTimestepUseHII == FALSE) ? -1 : dtPhoton;
  dtPhoton = 10*huge_number;

  FLOAT TimeNow = LevelArray[level]->GridData->ReturnTime();
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
	   &VelocityUnits, TimeNow);

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(TimeNow, &a, &dadt);
  float afloat = float(a);

  // Calculate timestep by limiting to a max change in HII
  if (RadiativeTransferHIIRestrictedTimestep) {

    for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
      for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
	ThisPhotonDT = Temp->GridData->
	  ComputePhotonTimestepHII(DensityUnits, LengthUnits, VelocityUnits, 
				   afloat, MetaData->GlobalMaximumkphIfront);
	dtPhoton = min(dtPhoton, ThisPhotonDT);
      }
    dtPhoton = CommunicationMinValue(dtPhoton);

    /* Use the average because the minimum ionization timescale can
       fluctuate significantly.  It gets even worse if the dtPhoton is
       allowed to vary a lot (>factor of a few). */
    
    //    printf("dtPhoton=%g, LastPhotonDT=%g \n", 
    //	   dtPhoton, LastPhotonDT); 

    if (LastPhotonDT > 0 && dtPhoton < huge_number) {// && LastTimestepUseHII) {
      if (dtPhoton > (1.0+MaxDTChange)*LastPhotonDT)
	dtPhoton = LastPhotonDT;
      else
	dtPhoton = 0.5 * (dtPhoton + LastPhotonDT);
    }

    if (dtPhoton < huge_number) {
      // Store dtPhoton before modifying it based on the next topgrid timestep
      LastPhotonDT = dtPhoton;  
      //LastTimestepUseHII = TRUE;
    }

  } // ENDIF

  // if we didn't find any cells that restrict timestep or the option
  // isn't requested, use hydro timestep on finest level
  if (dtPhoton >= huge_number) {
    if (maxLevel > 0) {
      for (Temp = LevelArray[maxLevel]; Temp; Temp = Temp->NextGridThisLevel) {
	ThisPhotonDT = Temp->GridData->ComputePhotonTimestep();
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
    //LastTimestepUseHII = FALSE;
  } // ENDIF

  /* Do not go past the level-0 time (+timestep) */

  HydroTime = LevelArray[0]->GridData->ReturnTime();
  if (level == 0)
    HydroTime += LevelArray[0]->GridData->ReturnTimeStep();
  float dtTol = PFLOAT_EPSILON * HydroTime;
  if ((HydroTime+dtTol - PhotonTime) < dtPhoton) {
    dtPhoton = max(HydroTime+dtTol - PhotonTime, dtTol);
    //LastTimestepUseHII = FALSE;
  }

  //LastTimestepUseHII = CommunicationMaxValue(LastTimestepUseHII);

  //if (InitialTimestep && !MetaData->FirstTimestepAfterRestart)
  //  dtPhoton = min(dtPhoton, dtLevelAbove);

  if (!RadiativeTransferAdaptiveTimestep && debug)
    printf("RadiativeTransfer: Setting dtPhoton = %g = %g years\n",
	   dtPhoton, dtPhoton*TimeUnits/3.1557e7);
  
  return SUCCESS;

}
