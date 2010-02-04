/***********************************************************************
/
/  REMOVE A PARTICLE BY ITS POSITION AND TYPE
/
/  written by: John Wise
/  date:       June, 2006
/  modified1:
/
/  RETURNS: 0 = not found; 1 = removed
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
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


int grid::RemoveParticle(int ID, bool disable)
{

  int i, found = FALSE;

  if (MyProcessorNumber != ProcessorNumber)
    return found;

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleNumber[i] == ID) {
      if (disable) {
	ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
	ParticleMass[i] = tiny_number;
      } else
	ParticleMass[i] = FLOAT_UNDEFINED;
      found = TRUE;
      break;
    }

  return found;

}
