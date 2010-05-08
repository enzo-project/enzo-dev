/***********************************************************************
/
/  GRID CLASS (COUNT STARS AND MIN LIFETIME)
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int grid::ReturnStarStatistics(int &Number, float &minLife)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int i, MassConvert;

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleType[i] != PARTICLE_TYPE_DARK_MATTER &&
	ParticleType[i] != PARTICLE_TYPE_TRACER &&
	ParticleType[i] != PARTICLE_TYPE_MUST_REFINE) {
      if (this->Time < ParticleAttribute[0][i]+ParticleAttribute[1][i])
	minLife = min(minLife, ParticleAttribute[1][i]);
      Number++;
    }

  return SUCCESS;

}
