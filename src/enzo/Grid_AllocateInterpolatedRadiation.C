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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::AllocateInterpolatedRadiation() 
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  int i,j,k, index, field, dim, vcsize = 1;

  for (dim = 0; dim < GridRank; dim++)
    vcsize *= GridEndIndex[dim] - GridStartIndex[dim] + 2;

  /* Initialize interpolated radiation fields if needed */

  int rkph, rgamma;
  
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (FieldsToInterpolate[field] == TRUE) {
      rkph = FieldUndefined;
      rgamma = FieldUndefined;
      switch (FieldType[field]) {
      case HIDensity:
	rkph = kphHINum;
	rgamma = gammaNum;
	break;
      case HeIDensity:
	if (RadiativeTransferHydrogenOnly == FALSE)
	  rkph = kphHeINum;
	break;
      case HeIIDensity:
	if (RadiativeTransferHydrogenOnly == FALSE)
	  rkph = kphHeIINum;
	break;
      } // ENDSWITCH field
      if (InterpolatedField[rkph] == NULL && rkph != FieldUndefined)
	InterpolatedField[rkph] = new float[vcsize];
      if (InterpolatedField[rgamma] == NULL && rgamma != FieldUndefined)
	InterpolatedField[rgamma] = new float[vcsize];
      for (i = 0; i < vcsize; i++)
	InterpolatedField[rkph][i] = 0.0f;
      if (rgamma != FieldUndefined)
	for (i = 0; i < vcsize; i++)
	  InterpolatedField[rgamma][i] = 0.0f;

    } // ENDIF fields to interpolate

  return SUCCESS;

}
