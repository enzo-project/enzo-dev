/***********************************************************************
/
/  GRID CLASS (INITIALIZE GRAVITATING MASS FIELD PARTICLES)
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
int grid::InitializeGravitatingMassFieldParticles(int RefinementFactor)
{
 
  /* Check to see if the field was already initialized. */
 
  if (GravitatingMassFieldParticlesCellSize != FLOAT_UNDEFINED)
    return SUCCESS;
 
  int dim, GravityBufferSize = GRAVITY_BUFFER_SIZE;
 
  /* Determine the size of the mass grid we'll need.  */
 
  for (dim = 0; dim < GridRank; dim++)
    switch (GravityBoundaryType) {
 
      /* TopGrid Periodic and TopGrid Isolated. */
 
    case TopGridPeriodic:
      GravityBufferSize = 0;
 
    case TopGridIsolated:
 
      /* Subgrid */
 
    case SubGridIsolated:
 
      /* Make the GravitatingMassField the size of the active region
	 plus the GravityBufferSize (in Parent cell units) on either size. */
 
      GravitatingMassFieldParticlesDimension[dim] =
	(GridEndIndex[dim] - GridStartIndex[dim] + 1) +
	  2*max(RefinementFactor*GravityBufferSize, NumberOfGhostZones);
      GravitatingMassFieldParticlesCellSize = CellWidth[dim][0];
      GravitatingMassFieldParticlesLeftEdge[dim] = GridLeftEdge[dim] -
	max(RefinementFactor*GravityBufferSize, NumberOfGhostZones)*
	  GravitatingMassFieldParticlesCellSize;
      break;
 
      /* Undefined or unknown is an error */
 
    case GravityUndefined:
    default:
      ENZO_FAIL("GravityBoundaryType undefined.\n");

    }
 
  /* Set unused dims. */
 
  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    GravitatingMassFieldParticlesDimension[dim] = 1;
 
  return SUCCESS;
}
