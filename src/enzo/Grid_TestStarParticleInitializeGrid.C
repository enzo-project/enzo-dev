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

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass,
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 int NumberOfTestStars, float clusterRadius)
{
  /* declarations */

  float CentralMass = 1.0;
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

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;

  printf("Central Mass: %f \n",CentralMass);

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = NumberOfTestStars;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }
  float p1;

  /* Set central particle. */
  for (i = 0; i <= NumberOfParticles; i++){
    for (dim = 0; dim < GridRank; dim++) {
      if (NumberOfParticles == 1){
        p1 = 0.5;
      }else{
        int rng = clusterRadius*200;
        p1 = float(rand() % rng)/100.0+(0.5-clusterRadius);
      }
        ParticlePosition[dim][i] = p1*
        (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
      ParticleVelocity[dim][i] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
    }
    ParticleMass[i] = CentralMass;
    ParticleAttribute[0][i] = Time+1e-7; //creation time:make sure it is non-zero
    if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][i] = 10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MOM_STAR))
      if(StarMakerExplosionDelayTime >= 0.0)
        ParticleAttribute[1][i] = 1.0;
      else
        ParticleAttribute[1][i] =10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MECHANICAL)) {
      if (StarParticleRadiativeFeedback){
        ParticleAttribute[1][i] = 25 * Myr_s/TimeUnits; // radiate for 25 Myr
      }
      else{
      ParticleAttribute[1][i] = 0.0;
      }
    }
  ParticleAttribute[2][i] = 0.0;  // Metal fraction
  ParticleAttribute[3][i] = 0.0;  // metalfSNIa
  }
  return SUCCESS;
}

