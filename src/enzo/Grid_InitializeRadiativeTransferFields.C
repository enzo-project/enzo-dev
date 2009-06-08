/***********************************************************************
/
/  GRID CLASS (INITIALIZE RADIATIVE TRANSFER FIELDS)
/
/  written by: Tom Abel
/  date:       June, 2003
/  modified1:
/
/ PURPOSE: we initialize photo and heating rates and the acceleration
/          field from radiation pressure.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::InitializeRadiativeTransferFields() 
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

//  if (HasRadiation == FALSE)
//    return SUCCESS;

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return FAIL;
  }

  int i,j,k, index;

  /* Initialize photo and heating rates and compute number densities */
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] =
	  BaryonField[gammaHINum][index]   =
	  BaryonField[kphHeINum][index]    = 
	  BaryonField[gammaHeINum][index]  =
	  BaryonField[kphHeIINum][index]   = 
	  BaryonField[gammaHeIINum][index] = 0.;
	if (MultiSpecies > 1 && !RadiativeTransferOpticallyThinH2) 
	  BaryonField[kdissH2INum][index] = 0.;
      }
    }  // loop over grid

  if (RadiationPressure) {

    int RPresNum1, RPresNum2, RPresNum3;
    if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) 
	== FAIL) {
      fprintf(stdout, "Error in IdentifyRadiationPressureFields.\n");
      return FAIL;
    }

    /* Initialize acceleration fields from radiation pressure */
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[RPresNum1][index] =
	    BaryonField[RPresNum2][index] =
	    BaryonField[RPresNum3][index] = 0.0;
	}
      }  // loop over grid
    
  }  /* ENDIF RadiationPressure */

  HasRadiation = FALSE;

  return SUCCESS;
}
