/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GRAVITY TEST)
/
/  written by: Greg Bryan
/  date:       April, 1995
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
 
 
int grid::TestGravityMotionInitializeGrid(float InitialVelocity)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int dim, i, NumberOfNewParticles = 1;
 
  /* Set particles. */
 
  if (NumberOfNewParticles > 0) {
 
    /* Set number of particles for this grid. */
 
    NumberOfParticles = NumberOfNewParticles;
 
    /* Allocate space. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    /* Set up particle positon & velocity. */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < NumberOfParticles; i++) {
	ParticlePosition[dim][i] = 0.5*(DomainLeftEdge[dim] +
					DomainRightEdge[dim]);
	ParticleVelocity[dim][i] = 0;
	if (dim == 0)
	  ParticleVelocity[dim][i] = InitialVelocity;
	ParticleMass[i]          = 1.0;
	ParticleNumber[i]        = i;
	ParticleType[i]          = PARTICLE_TYPE_DARK_MATTER;
      }
 
  }
 
  return SUCCESS;
}
