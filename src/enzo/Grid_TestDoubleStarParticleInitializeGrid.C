/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestDoubleStarParticleInitializeGrid(FLOAT TestStarParticleStarMass[2], 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[2][3],
					 FLOAT TestStarParticleStarPosition[2][3],
					 FLOAT TestStarParticleMetallicity[2])
{
  /* declarations */

  float StarMass1 = 1.0;
  float StarMass2 = 1.0;
  int i, dim;
  float TestInitialdt = *Initialdt;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = 2;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }

  /* Set central particle. */ 
  for (i = 0; i < NumberOfParticles; i++) {
    for (dim = 0; dim < GridRank; dim++) {
        ParticlePosition[dim][i] = TestStarParticleStarPosition[i][dim];
        ParticleVelocity[dim][i] = TestStarParticleStarVelocity[i][dim]*1e5*TimeUnits/LengthUnits;
    }
    ParticleMass[i] = TestStarParticleStarMass[i]*1.99e33* POW(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;
    if (StarMakerStoreInitialMass)
      ParticleInitialMass[i] = TestStarParticleStarMass[i]*1.99e33* POW(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;
    ParticleAttribute[0][i] = Time+1e-7;

    ParticleAttribute[1][i] = StarMakerMinimumDynamicalTime * yr_s/TimeUnits;
    if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][i] = 10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MOM_STAR) || STARFEED_METHOD(MECH_STAR))
        if(StarMakerExplosionDelayTime >= 0.0)
        ParticleAttribute[1][i] = 1.0;
    
    ParticleAttribute[2][i] = TestStarParticleMetallicity[i];  // Metal fraction
    ParticleAttribute[3][i] = 0.0;  // metalfSNIa
  }

  return SUCCESS;
}

