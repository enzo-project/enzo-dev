/***********************************************************************
/
/  GRID CLASS (CONVERTS DARK MATTER PARTICLES TO MUST REFINE IN REGION)
/
/  written by: Christine Simpson & Greg Bryan
/  date:       October, 2009
/  modified1:
/
/  PURPOSE:Flag particles within a rectangular region as must refine 
/  particles.  The region's dimensions are set by the parameters 
/  MustRefineRegionLeftEdge and MustRefineRegionRightEdge.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::MustRefineParticlesFlagInRegion()
{
  /* Local variable declarations. */

  int n, dim, INSIDE_REGION, NumberOfParticlesConverted = 0;

  /* Exit if grid is not on this processor, but not before recording increase
     in the number of particles, both in this grid and in the metadata. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Go through all the particles, and if they are inside the MustRefineRegion,
     convert from DARK_MATTER_PARTICLE to MUST_REFINE_PARTICLE. */

  for (n = 0; n < NumberOfParticles; n++) {
    if (ParticleType[n] == PARTICLE_TYPE_DARK_MATTER) {
      INSIDE_REGION = TRUE;
      for (dim = 0; dim < GridRank; dim++) {
	if (ParticlePosition[dim][n] < MustRefineParticlesLeftEdge[dim] ||
	    ParticlePosition[dim][n] > MustRefineParticlesRightEdge[dim]) {
	  INSIDE_REGION = FALSE;
	  break;
	}
      }
      if (INSIDE_REGION == TRUE) {
	ParticleType[n] = PARTICLE_TYPE_MUST_REFINE;	  
	NumberOfParticlesConverted++;
      }
    }
  }
  if (NumberOfParticlesConverted > 0)
    printf("MustRefineParticlesFlagInRegion:ParticlesFlagged = %d\n", NumberOfParticlesConverted);

  return SUCCESS;
}
