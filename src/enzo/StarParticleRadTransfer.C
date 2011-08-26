/***********************************************************************
/
/  CONVERT SHINE PARTICLE INTO RADIATIVE TRANSFER PARTICLE
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/ PURPOSE: This routine converts particles that shone by a 1/r^2 law
/          into particles that utilize an adaptive 3D ray tracing
/          scheme.
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

#define NO_MEAN_ENERGY

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleRadTransfer(LevelHierarchyEntry *LevelArray[], int level,
			    Star *AllStars)
{

  /* If photon test simulation, don't change the radiation sources. */

  if (ProblemType == 50)
    return SUCCESS;

  int i, j, nShine, nbins;
  double Q[MAX_ENERGY_BINS], QTotal;
  float qFrac, lramp, energies[MAX_ENERGY_BINS];
  float XRayLuminosityFraction = 0.43;
  Star *cstar;

  /* If sources exist, delete them */

  RadiationSourceEntry *dummy;
  while (GlobalRadiationSources != NULL) {
    dummy = GlobalRadiationSources;
    GlobalRadiationSources = GlobalRadiationSources->NextSource;
    delete dummy;
  }

  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;

  if (AllStars == NULL)
    return SUCCESS;

  /* Retrieve the units */

  FLOAT Time = LevelArray[level]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);

  // Convert to years
  float TimeInYears = 3.1557e7 / TimeUnits;

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

    // Check the rules if this star particle is radiative
    if (cstar->IsARadiationSource(Time)) {

      // Calculate photon luminosity
      if (cstar->ComputePhotonRates(nbins, energies, Q) == FAIL) {
	ENZO_FAIL("Error in ComputePhotonRates.\n");
      }
      
      QTotal = 0;
      for (j = 0; j < nbins; j++) QTotal += Q[j];
      for (j = 0; j < nbins; j++) Q[j] /= QTotal;

#ifdef USE_MEAN_ENERGY
      double meanEnergy = 0;
      nbins = 1;
      for (j = 0; j < nbins; j++)
	meanEnergy += energies[j] * Q[j];
      meanEnergy /= QTotal;
      energies[0] = meanEnergy;
      Q[0] = QTotal;
#endif /* USE_MEAN_ENERGY */

      /* (TODO) If requested, calculate ramping time for the luminosity */

      float ramptime = 0.0;   // zero for no ramp
      float tdyn, ti;
      if (cstar->ReturnType() == PopII) {
	if (StarClusterUnresolvedModel) {  // Cen & Ostriker
	  ramptime = cstar->ReturnLifetime();
	} else {  // Wise & Cen
	  ramptime = TimeInYears * StarClusterMinDynamicalTime;
	}
      } else if (cstar->ReturnType() == PopIII)
	// should be an parameter or determined from the data
	ramptime = TimeInYears * 50e3;
      else if (cstar->ReturnType() == SimpleSource)
	ramptime = TimeInYears * 1e6 * SimpleRampTime;

      /* Transfer the shining particle properties to the radiative
	 transfer source particle */

      RadiationSourceEntry *RadSource;
      RadSource = cstar->RadiationSourceInitialize();
      RadSource->Luminosity     = QTotal * LConv;
      RadSource->RampTime       = ramptime;
      RadSource->EnergyBins     = nbins;
      RadSource->Energy         = new float[nbins];
      RadSource->SED            = new float[nbins];
      
      for (j = 0; j < nbins; j++) {
	RadSource->Energy[j] = energies[j];
	RadSource->SED[j]    = Q[j];
      }

      // if the source needs a beaming direction, define it here
      RadSource->Orientation    = NULL;

      if (GlobalRadiationSources->NextSource != NULL)

	GlobalRadiationSources->NextSource->PreviousSource = RadSource;
      GlobalRadiationSources->NextSource = RadSource;
      
    } // ENDIF is a radiation source?

  } // ENDFOR stars

  return SUCCESS;

}
