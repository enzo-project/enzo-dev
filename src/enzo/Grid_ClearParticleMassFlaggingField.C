/***********************************************************************
/
/  GRID CLASS (CLEAR THE PARTICLE MASS FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  May, 2009 by John Wise: for particles
/
/  PURPOSE:
/
************************************************************************/
 
// Allocate and clear the particle mass flagging field.
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::ClearParticleMassFlaggingField()
{
 
  /* error check */
 
  if (ParticleMassFlaggingField != NULL) {
    fprintf(stderr, "ClearParticleMassFlaggingField: Warning, field not deleted.\n");
    delete [] ParticleMassFlaggingField;
  }
 
  /* compute size and allocate */
 
  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  ParticleMassFlaggingField = new float[size];
 
  /* Clear it */
 
  for (i = 0; i < size; i++)
    ParticleMassFlaggingField[i] = 0.0;
 
}
