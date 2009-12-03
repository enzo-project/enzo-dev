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

#define RT_ENERGY_BINS 4
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

  int i, j, nShine;
  double Q[RT_ENERGY_BINS], QTotal;
  float qFrac, lramp, energies[RT_ENERGY_BINS];
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

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

    // Check the rules if this star particle is radiative
    if (cstar->IsARadiationSource(Time)) {

      // Calculate photon luminosity
      if (cstar->ComputePhotonRates(energies, Q) == FAIL) {
	fprintf(stderr, "Error in ComputePhotonRates.\n");
	ENZO_FAIL("");
      }

      QTotal = 0;
      for (j = 0; j < RT_ENERGY_BINS; j++) QTotal += Q[j];
      for (j = 0; j < RT_ENERGY_BINS; j++) Q[j] /= QTotal;

      int nbins = RT_ENERGY_BINS;
      double meanEnergy = 0;

#ifdef USE_MEAN_ENERGY
      nbins = 1;
      for (j = 0; j < RT_ENERGY_BINS; j++)
	meanEnergy += energies[j] * Q[j];
      meanEnergy /= QTotal;
      energies[0] = meanEnergy;
      Q[0] = QTotal;
#endif /* USE_MEAN_ENERGY */

      /* (TODO) If requested, calculate ramping time for the luminosity */

      float ramptime = 0.0;   // zero for no ramp

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

      GlobalRadiationSources->NextSource = RadSource;
      
    } // ENDIF is a radiation source?

  } // ENDFOR stars

  return SUCCESS;

}
