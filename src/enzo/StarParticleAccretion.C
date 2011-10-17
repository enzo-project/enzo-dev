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
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdlib.h>
#include <stdio.h>
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
#include "CommunicationUtilities.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RecalibrateAccretingMass(FLOAT star_pos[], LevelHierarchyEntry *LevelArray[], 
			     int level, float &BondiRadius, float &density,
			     float &RecalibrateAccretingMassRatio);

int StarParticleAccretion(TopGridData *MetaData, 
			  LevelHierarchyEntry *LevelArray[], int level, 
			  Star *&AllStars)
{

#define NOT_SEDOV_TEST
#define NOT_HII_REGION_TEST

#if defined(SEDOV_TEST) || defined(HII_REGION_TEST)
  return SUCCESS;
#endif

  LCAPERF_START("StarParticleAccretion");

  Star *ThisStar;
  FLOAT Time;
  float BondiRadius = 0.0, density = 0.0, RecalibrateAccretingMassRatio = 1.0;
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

    // needed for RecalibrateAccretingMass
    BondiRadius = 0.0; density = 0.0;

    /* Calculate accreting mass */

    if (ThisStar->CalculateMassAccretion(BondiRadius, density) == FAIL) {
      fprintf(stderr, "Error in star::CalculateMassAccretion.\n");
      ENZO_FAIL("");
    }
    
    /* Correct the accreting mass if requested */

    if ((MBHAccretion > 0) &&
	(MBHAccretingMassRatio == BONDI_ACCRETION_CORRECT_NUMERICAL)) {

      // communicate BondiRadius and density as it is known on only one processor
      CommunicationAllSumValues(&BondiRadius, 1);
      CommunicationAllSumValues(&density, 1);

      // calibrate the accretion rate
      RecalibrateAccretingMassRatio = 1.0;
      RecalibrateAccretingMass(ThisStar->ReturnPosition(), LevelArray, level, 
      			       BondiRadius, density, RecalibrateAccretingMassRatio);
      fprintf(stdout, "BondiRadius = %g, RecalibrateAccretingMassRatio = %g\n", 
	      BondiRadius, RecalibrateAccretingMassRatio);

      if ((BondiRadius > 0.0) && (ThisStar->ReturnType() == MBH) && 
	  (ThisStar->ReturnCurrentGrid() != NULL)) 
	ThisStar->MultiplyAccretionRate(RecalibrateAccretingMassRatio);

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

  LCAPERF_STOP("StarParticleAccretion");
  return SUCCESS;

}
