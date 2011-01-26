/***********************************************************************
/
/  GRID CLASS (ALLOCATE AND CLEAR THE GRAVITATING MASS FIELD PARTICLES)
/
/  written by: Greg Bryan
/  date:       June, 1995
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
 
 
int grid::ClearGravitatingMassFieldParticles()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (GravitatingMassFieldParticlesCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("GravitatingMassFieldParticles uninitialized.\n");
  }
 
  /* Compute size of the gravitating mass field. */
 
  int dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];
 
  /* Allocate and clear the field. */
 
//  if (GravitatingMassFieldParticles != NULL)
//    fprintf(stderr, "ClearGravitatingMassField: Warning! Field not NULL.\n");
 
  if (GravitatingMassFieldParticles == NULL)
    GravitatingMassFieldParticles = new float[size];
  if (GravitatingMassFieldParticles == NULL) {
    ENZO_FAIL("malloc error (out of memory?)\n");

  }
 
  for (int i = 0; i < size; i++)
    GravitatingMassFieldParticles[i] = 0.0;
 
  return SUCCESS;
}
