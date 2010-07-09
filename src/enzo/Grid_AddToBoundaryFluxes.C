/***********************************************************************
/
/  GRID CLASS (ADD FLUXES IN ARGUMENT TO THE BOUNDARY FLUXES OF THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Solve the hydro equations with the solver, saving the subgrid fluxes
//
 
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
 
int grid::AddToBoundaryFluxes(fluxes *BoundaryFluxesToBeAdded)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro)
    return SUCCESS;
 
  int i, j, dim, field;
 
  /* Error check */
 
  if (BoundaryFluxes == NULL) {
    ENZO_FAIL("grid->AddToBoundarFluxes BoundaryFluxes not defined.\n");
  }

  if (BoundaryFluxesToBeAdded == NULL) {
    ENZO_FAIL("grid->AddToBoundarFluxes BoundaryFluxesToBeAdded not defined.\n");

  }

 
  /* add the fluxes of the structure pointed to by the argument to
     the fluxes of this grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* compute size (in floats) of flux storage and then copy.  */
 
    int size = 1;
    for (j = 0; j < GridRank; j++)
      size *=
	BoundaryFluxesToBeAdded->LeftFluxEndGlobalIndex[dim][j] -
	  BoundaryFluxesToBeAdded->LeftFluxStartGlobalIndex[dim][j] + 1;
 
    /* copy */
 
    for (field = 0; field < NumberOfBaryonFields; field++)
      for (i = 0; i < size; i++) {
	BoundaryFluxes->LeftFluxes[field][dim][i] +=
	  BoundaryFluxesToBeAdded->LeftFluxes[field][dim][i];
	BoundaryFluxes->RightFluxes[field][dim][i] +=
	  BoundaryFluxesToBeAdded->RightFluxes[field][dim][i];
      }
  }
 
  return SUCCESS;
 
}
