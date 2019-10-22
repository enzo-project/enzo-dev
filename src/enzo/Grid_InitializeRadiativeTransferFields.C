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
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::InitializeRadiativeTransferFields() 
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

//  if (HasRadiation == FALSE)
//    return SUCCESS;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kdissH2IINum, kphHMNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  int PeNum, FUVRateNum;

  PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);
  FUVRateNum = FindField(FUVRate, this->FieldType, this->NumberOfBaryonFields);


  int RaySegNum = FindField(RaySegments, FieldType, NumberOfBaryonFields);

  int i,j,k, index;

  /* Initialize photo and heating rates and compute number densities */
  /* All the fields should be zeroed - even if all won't be used.... */
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] = BaryonField[gammaNum][index]   = 0.0;
      }
      
    }  // loop over grid

  if (RadiativeTransferHydrogenOnly == FALSE) 
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kphHeINum][index]  = BaryonField[kphHeIINum][index] = 0.0;
	}
      }  // loop over grid

  if (MultiSpecies > 1 && !RadiativeTransferFLD)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kdissH2INum][index] = BaryonField[kphHMNum][index] = 0.0;
	  BaryonField[kdissH2IINum][index] = 0.0;
	}
      }  // loop over grid

  if (RadiationPressure) {

    int RPresNum1, RPresNum2, RPresNum3;
    IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3);

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

  if (RadiativeTransferLoadBalance)
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = 0; i < GridDimension[0]; i++, index++)
	  BaryonField[RaySegNum][index] = 0.0;
      }  // loop over grid

  /* This is a litte rough for now, but going through testing -
     this can be distilled to just a couple params once
     done - Dec 2016 - AJE
  */

  TIMER_START("InitializeOTFields");
  if (RadiativeTransferOpticallyThinFUV    &&
      IndividualStarFUVHeating){
    if (PeNum < 0) ENZO_FAIL("Failure to identify PeHeatingRate in InitializeRadiativeTransferFields\n");
    for ( k = 0; k < GridDimension[2]; k++){
      for( j = 0; j < GridDimension[1]; j++){
        index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, index++){
          BaryonField[PeNum][index] = 0.0;
        }
      }
    }
  }

  if (!RadiativeTransferOpticallyThinFUV   &&
      IndividualStarFUVHeating){
    for ( k = 0; k < GridDimension[2]; k++){
      for( j = 0; j < GridDimension[1]; j++){
        index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, index++){
          BaryonField[PeNum][index] = 0.0;
          BaryonField[FUVRateNum][index] = 0.0;
        }
      }
    }
  }

// -- end AJE comments

  TIMER_STOP("InitializeOTFields");

  HasRadiation = FALSE;
  MaximumkphIfront = 0;

  return SUCCESS;
}
