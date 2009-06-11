/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (CONSTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
//
// Constructor
//
 
ExternalBoundary::ExternalBoundary()
{
 
  /* declarations */
 
  int dim, field, i;
 
  /* Set scalars to zero. */
 
  BoundaryRank         = 0;
  NumberOfBaryonFields = 0;
  ParticleBoundaryType = BoundaryUndefined;
 
  /* Set all dims to 1. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    BoundaryDimension[dim] = 1;
 
  /* Clear BoundaryType and BoundaryValue pointers. */
 
  for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++)
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      for (i = 0; i < 2; i++) {
	BoundaryType[field][dim][i] = NULL;
	BoundaryValue[field][dim][i] = NULL;
      }
	
}
