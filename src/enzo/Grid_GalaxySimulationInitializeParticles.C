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


  // remove this if need to sep by processor
//   NumberOfParticles = NumberOfDMParticles;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

//  if (this->NumberOfSubgrids > 1) return SUCCESS;

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

  const float msolar = 1.989E33;

  int count = 0;
  FLOAT cell_volume = this->CellWidth[0][0]*this->CellWidth[0][0]*this->CellWidth[0][0];

  int *particle_index;
  particle_index = new int[NumberOfDMParticles];

  /* find which particles belong in the simulation
     handle unit and coordinate conversions */
  for(int i = 0; i < NumberOfDMParticles; i++){

    int off_grid = 0;
    for(int dim = 0; dim < MAX_DIMENSION; dim ++){
      DMParticlePosition[dim][i] = DMParticlePosition[dim][i]*pc/LengthUnits + DiskGravityPosition[dim];

      off_grid += !( (DMParticlePosition[dim][i] > this->CellLeftEdge[dim][NumberOfGhostZones]   )*
                     (DMParticlePosition[dim][i] < this->CellLeftEdge[dim][this->GridDimension[dim] - NumberOfGhostZones]) );
    }

    particle_index[i] = -1;

   // off_grid = FALSE; // remove if need to separate by processor

    if (off_grid) continue;
    particle_index[count] = i;
    count++;
  }

  /* allocate memory for the new particles */
  if (count > 0 )
    this->AllocateNewParticles(count);


  /* now deposit the particles that belong to this grid */
  if (count > 0){
    for(int i = 0; i < count; i ++){
      for (int dim = 0; dim < MAX_DIMENSION; dim ++){
        ParticlePosition[dim][i] = DMParticlePosition[dim][particle_index[i]];
        ParticleVelocity[dim][i] = DMParticleVelocity[dim][particle_index[i]] * 1.0E5 / VelocityUnits;
      }
   
      ParticleMass[i] = DMParticleMass[particle_index[i]] * msolar / MassUnits / (cell_volume);

      ParticleAttribute[0][i] = this->Time;
      ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
      ParticleNumber[i] = particle_index[i]; // set in GalaxySimulationInitialize
    }
  }

  printf("P(%"ISYM"): Deposited %"ISYM" Dark Matter Particles out of %"ISYM"\n", MyProcessorNumber, count, NumberOfDMParticles);

  if ((NumberOfDMParticles - count) > 0){
    printf("P(%"ISYM"): WARNING: %"ISYM" particles were not within the computational domain \n", MyProcessorNumber, NumberOfDMParticles - count);
  }

  this->NumberOfParticles = count;

  // clean up
  delete [] particle_index;

  return SUCCESS;
}
