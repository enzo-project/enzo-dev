/***********************************************************************
/
/  DELETE FLUXES FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// Delete all the data associated with the fluxes structure in the argument
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "Fluxes.h"
 
void DeleteFluxes(fluxes *Fluxes)
{
  if (Fluxes != NULL)
    for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++)
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	if (Fluxes->LeftFluxes[field][dim] != NULL)
	  delete [] Fluxes->LeftFluxes[field][dim];
	if (Fluxes->RightFluxes[field][dim] != NULL)
	  delete [] Fluxes->RightFluxes[field][dim];
	Fluxes->LeftFluxes[field][dim]  = NULL;
	Fluxes->RightFluxes[field][dim] = NULL;
      }
}
