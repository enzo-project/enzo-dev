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
#include "StarParticleData.h"

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

int grid::RemoveParticle(int ID)
{

  int i, found = FALSE;

  if (MyProcessorNumber != ProcessorNumber)
    return found;

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleNumber[i] == ID) {
      //      ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
      //      ParticleMass[i] = tiny_number;
      ParticleMass[i] = FLOAT_UNDEFINED;
      found = TRUE;
      break;
    }

  return found;

}
