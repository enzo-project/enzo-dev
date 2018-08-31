/***********************************************************************
/
/  GRID CLASS (COMPUTE CENTER OF MASS OF ALL PARTICLES ON THIS GRID)
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

int grid::DiskGravityComputeParticleCOM(FLOAT *localCOM, float & localMass){

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  FLOAT x = 0.0, y = 0.0, z = 0.0;
  float m = 0.0, weight = 0.0;

  for(int i = 0; i < this->NumberOfParticles; i++){
    weight = this->ParticleMass[i];

    x     += weight * this->ParticlePosition[0][i];
    y     += weight * this->ParticlePosition[1][i];
    z     += weight * this->ParticlePosition[2][i];
    m     += weight;
  }

  // now compute grid COM and mass on this grid
  float inv_m = 0.0;
  if (m > 0) 
    inv_m = 1.0 / m;

  localCOM[0] = x * inv_m;
  localCOM[1] = y * inv_m;
  localCOM[2] = z * inv_m;
  localMass   = m;

  return SUCCESS;
}
