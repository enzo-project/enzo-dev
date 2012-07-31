/***********************************************************************
/
/  GRID CLASS (SEARCH FOR ALL STAR PARTICLES AND RETURN HOW MANY)
/
/  written by: Greg Bryan
/  date:       September, 2000
/  modified1:  JHK & JHW (2009)
/
/  PURPOSE:
/
/  NOTE:
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

int grid::ReturnNumberOfStarParticles(void)
{
  int np = 0;
  int n, abstype;
  if (MyProcessorNumber == ProcessorNumber)
    for (n = 0; n < NumberOfParticles; n++) {
      abstype = ABS(ParticleType[n]);
      if (abstype == PARTICLE_TYPE_STAR ||
	  (abstype >= PARTICLE_TYPE_MUST_REFINE &&
	   abstype != PARTICLE_TYPE_MBH)) np++;
    }
  return np;
}
