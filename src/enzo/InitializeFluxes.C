/***********************************************************************
/
/  INITIALIZE FLUXES FUNCTION
/
/  written by: Sam Skillman 
/  date:       October, 2013 
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// Initialize all the data associated with the fluxes structure in the argument
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "Fluxes.h"
 
void InitializeFluxes(fluxes *Fluxes)
{

  for (int i=0; i<MAX_DIMENSION; i++){
    for (int j=0; j<MAX_DIMENSION; j++){
      Fluxes->LeftFluxStartGlobalIndex[i][j] = 0;
      Fluxes->LeftFluxEndGlobalIndex[i][j] = 0;
      Fluxes->RightFluxStartGlobalIndex[i][j] = 0;
      Fluxes->RightFluxEndGlobalIndex[i][j] = 0;
    }
  }

  for (int i=0; i<MAX_NUMBER_OF_BARYON_FIELDS; i++){
    for (int j=0; j<MAX_DIMENSION; j++){
      Fluxes->LeftFluxes[i][j] = NULL;
      Fluxes->RightFluxes[i][j] = NULL;
    }
  }
}
