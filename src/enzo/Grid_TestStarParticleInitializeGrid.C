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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[])
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

  NumberOfParticles = 1;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }

  /* Set central particle. */ 
  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = TestStarParticleStarPosition[dim]*
      (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
    ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
  }
  ParticleMass[0] = CentralMass;
  ParticleAttribute[0][0] = Time+1e-7; //creation time:make sure it is non-zero
  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 1e7*3.15e7/TimeUnits;
  if (STARFEED_METHOD(MOM_STAR)) ParticleAttribute[1][0] = TestInitialdt/100.0; // dynamical time?
  ParticleAttribute[2][0] = 0.0;  // Metal fraction
  ParticleAttribute[3][0] = 0.0;  // metalfSNIa

  return SUCCESS;
}

