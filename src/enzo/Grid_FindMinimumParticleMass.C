/***********************************************************************
/
/  FINDS THE MINIMUM PARTICLE MASS
/
/  written by: John Wise
/  date:       April, 2009
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

int grid::FindMinimumParticleMass(float &min_mass, int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i;
  float MassFactor = pow(8.0, -level);

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleType[i] == PARTICLE_TYPE_DARK_MATTER)
      min_mass = min(min_mass, MassFactor*ParticleMass[i]);

  return SUCCESS;

}
