/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A SELFFORCE TEST)
/
/  written by: Jean-Claude Passy
/  date:       June 2013
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
 
 
int grid::TestSelfForceInitializeGrid(float CentralDensity,
				      int NumberOfNewParticles,
				      float xpos, float ypos, float zpos,
				      float vx, float vy, float vz)
{
  /* declarations */
 
  int dim, i, size, field, vel;
  float phi, r, theta, pi = 3.14159;
  
  /* Return if this doesn't concern us. */
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfNewParticles > 0) {

    /* Set particles. */
    NumberOfParticles = NumberOfNewParticles;
 
    /* Allocate space. */
    this->AllocateNewParticles(NumberOfParticles);
 
    /* Set central particle */
    ParticleMass[0] = CentralDensity;

    ParticlePosition[0][0] = xpos;
    ParticlePosition[1][0] = ypos;
    ParticlePosition[2][0] = zpos;

    ParticleVelocity[0][0] = vx;
    ParticleVelocity[1][0] = vy;
    ParticleVelocity[2][0] = vz;
  }

  return SUCCESS;
}
