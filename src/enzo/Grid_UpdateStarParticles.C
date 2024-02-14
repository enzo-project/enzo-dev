/***********************************************************************
/
/  CREATES STAR PARTICLES FROM EXISTING PARTICLES
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:  negative types mark particles that have just been before 
/          and not been converted into a star particle.
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

int grid::UpdateStarParticles(int level, std::map<int, Star*>* const &StarParticleLookupMap)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0 || Stars == NULL)
    return SUCCESS;

  int i;
  Star *cstar;

// this bit of code uses a lookup map to find and copy star particles,
// saving computation when there are many star particles
  for (i = 0; i < NumberOfParticles; i++) {
    auto search = (*StarParticleLookupMap).find(ParticleNumber[i]);
    if (search != (*StarParticleLookupMap).end()) {
      cstar = search->second;
        if (cstar->type > 0) // living stars only (<0 == waiting to be created)
          cstar->CopyFromParticle(this, i, level);
    }
  }

  return SUCCESS;

}
