/***********************************************************************
/
/  REFLECTS CHANGES IN STAR PARTICLE IN NORMAL PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
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
#include "Star.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::MirrorStarParticles(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  const double Msun = 1.989e33;
  int i, dim, place;
  float MassConversion;
  Star *cstar;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  MassConversion = (float) (double(LengthUnits*CellWidth[0][0]) *
			    double(LengthUnits*CellWidth[0][0]) *
			    double(LengthUnits*CellWidth[0][0]) *
			    double(DensityUnits) / Msun);

  for (cstar = Stars; cstar; cstar = cstar->NextStar) {

    // Find where this star particle is stored in main arrays
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleNumber[i] == cstar->Identifier) {
	place = i;
	break;
      }

    // Change all particle data in favor of updated Star particle
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ParticlePosition[dim][place] = cstar->pos[dim];
      ParticleVelocity[dim][place] = cstar->vel[dim];
    }
    ParticleMass[place] = (float) (cstar->Mass / MassConversion);
    ParticleType[place] = cstar->type;
    ParticleAttribute[0][place] = cstar->BirthTime;
    ParticleAttribute[1][place] = cstar->LifeTime;

  } // ENDFOR stars
  
  return SUCCESS;
}
