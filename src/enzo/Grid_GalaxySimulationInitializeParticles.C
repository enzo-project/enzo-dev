/***********************************************************************
/
/  GRID CLASS (INITIALIZE DARK MATTER PARTICLES IN GALAXY SIMULATION)
/
/  written by: Andrew Emerick
/  date:       March, 2017
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
             float *VelocityUnits, FLOAT Time);

int grid::GalaxySimulationInitializeParticles(int NumberOfDMParticles, 
                                              float *DMParticleMass, FLOAT *DMParticlePosition[],
                                              float *DMParticleVelocity[]){

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (this->NumberOfSubgrids > 1) return SUCCESS;

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

  const float msolar = 1.989E33;

  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  int count = 0;
  FLOAT cell_volume = this->CellWidth[0][0]*this->CellWidth[0][0]*this->CellWidth[0][0];

  /* now loop through and deposit particles if they are on this grid */
  for(int i = 0; i < NumberOfDMParticles; i++){

    int off_grid = 0;
    for(int dim = 0; dim < MAX_DIMENSION; dim ++){
      DMParticlePosition[dim][i] = DMParticlePosition[dim][i]*pc/LengthUnits + DiskGravityPosition[dim];

      off_grid += !( (DMParticlePosition[dim][i] > this->CellLeftEdge[0][ibuff])*
                     (DMParticlePosition[dim][i] < this->CellLeftEdge[dim][nx-ibuff]) );
    }

    if (off_grid) continue;

    for (int dim = 0; dim < MAX_DIMENSION; dim ++){
      ParticlePosition[dim][count] = DMParticlePosition[dim][i];
      ParticleVelocity[dim][count] = DMParticleVelocity[dim][i] * 1.0E5 / VelocityUnits;
    }
    ParticleMass[count] = DMParticleMass[i] * msolar / MassUnits / (cell_volume);

    ParticleAttribute[0][count] = this->Time;
    ParticleType[count] = PARTICLE_TYPE_DARK_MATTER;
    ParticleNumber[count] = count;

  }

  printf("P(%"ISYM"): Deposited %"ISYM" Dark Matter Particles\n", MyProcessorNumber, count);

  return SUCCESS;
}
