/***********************************************************************
/
/  ADD 1/R^2 RADIATION FIELD FOR SOURCES CREATED SINCE THE LAST FS 
/  SOLVER CALL
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/ PURPOSE: 
/          
/
************************************************************************/
#include "preincludes.h"

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
#include "StarParticleData.h"
#include "RadiativeTransferHealpixRoutines.h"
#include "ImplicitProblemABC.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
float CommunicationMinValue(float Value);

int FLDCorrectForImpulses(int field, LevelHierarchyEntry *LevelArray[],
			  int level, Star *AllStars, FLOAT FLDTime, 
			  float dtFLD)
{

  RadiationSourceEntry *RS = NULL;
  Star *cstar = NULL;

  int nbins;
  FLOAT TimeNow, BirthTime, *Position;
  double Luminosity, LConv, LL[MAX_ENERGY_BINS], sigma;
  float TimeFraction, Lifetime, energies[MAX_ENERGY_BINS];
  LevelHierarchyEntry *Temp;

  TimeNow = LevelArray[level]->GridData->ReturnTime();

  switch (field) {
  case kdissH2I:
    sigma = 3.71e-18;
    break;
  default:
    if (debug) fprintf(stderr, "radiation field = %"ISYM"\n", field);
    ENZO_FAIL("Free-streaming radiation field not recognized.");
    break;
  }

  /* Loop over stars or sources, looking for ones that were created
     after the last FLD call */

  if (ProblemType == 50)
    RS = GlobalRadiationSources->NextSource;
  else
    cstar = AllStars;

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Convert from #/s to RT units
  LConv = (double) TimeUnits / pow(LengthUnits,3);

  while (RS != NULL && cstar != NULL) {

    if (ProblemType == 50) {
      // Check if this source was created since the last FLD timestep.
      // Only add its flux once!
      if (RS->AddedEmissivity || RS->CreationTime < FLDTime-dtFLD)
	continue;
      BirthTime = RS->CreationTime;
      Position = RS->Position;
      Luminosity = RS->Luminosity * RS->SED[3] / LConv;
      TimeFraction = (min(RS->CreationTime + RS->LifeTime, FLDTime) -
		      max(RS->CreationTime, FLDTime-dtFLD)) / dtFLD;
      RS->AddedEmissivity = true;
    } else {
      BirthTime = cstar->ReturnBirthTime();
      if (cstar->ReturnEmissivityFlag() || BirthTime < FLDTime-dtFLD)
	continue;
      Lifetime = cstar->ReturnLifetime();
      Position = cstar->ReturnPosition();
      cstar->ComputePhotonRates(nbins, energies, LL);
      Luminosity = LL[3];
      TimeFraction = (min(BirthTime + Lifetime, FLDTime) -
		      max(BirthTime, FLDTime-dtFLD)) / dtFLD;
      cstar->AddEmissivityFlag();
    }

    Luminosity *= TimeFraction;
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->AddRadiationImpulse(field, Luminosity, sigma, BirthTime,
					    Position);
    
    if (ProblemType == 50)
      RS = RS->NextSource;
    else
      cstar = cstar->NextStar;

  } // ENDWHILE sources/stars

  return SUCCESS;

}
