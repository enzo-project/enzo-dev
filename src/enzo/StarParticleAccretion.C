/***********************************************************************
/
/  ADDS ACCRETED MASS TO PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             September, 2009
/
/  NOTES:
/
************************************************************************/

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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleAccretion(TopGridData *MetaData, 
			  LevelHierarchyEntry *LevelArray[], int level, 
			  Star *&AllStars)
{

#define NOT_SEDOV_TEST
#define NOT_HII_REGION_TEST

#if defined(SEDOV_TEST) || defined(HII_REGION_TEST)
  return SUCCESS;
#endif

  Star *ThisStar;
  FLOAT Time;
  LevelHierarchyEntry *Temp;

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);


  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {

    // Must be the finest level for any feedback except star formation
    if (ThisStar->ReturnFeedbackFlag() != FORMATION &&
	LevelArray[level+1] != NULL)
      continue;

    if (ThisStar->CalculateMassAccretion() == FAIL) {
      ENZO_FAIL("Error in star::CalculateMassAccretion.\n");
    }

    /* Add accreted mass to star particles */

    if (ThisStar->Accrete() == FAIL) {
      ENZO_FAIL("Error in star::Accrete.\n");
    }

    /* Add accreted angular momentum to star particles */

    if (ThisStar->AccreteAngularMomentum() == FAIL) {
      ENZO_FAIL("Error in star::AccreteAngularMomentum.\n");

    }

  }

  return SUCCESS;

}
