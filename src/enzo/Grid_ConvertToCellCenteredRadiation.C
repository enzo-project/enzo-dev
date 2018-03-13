#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CELL CENTERED RADIATION FIELD FROM VERTEX-CENTERED)
/
/  written by: John Wise
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: Takes vertex-centered field in InterpolatedField[] and 
/           computes the cell-centered counterpart are puts it in 
/           BaryonField[].  Must divide by cell-centered absorbed 
/           density.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"

int grid::ConvertToCellCenteredRadiation()
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (!HasRadiation)
    return SUCCESS;

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  int rkph, rgamma, field;

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

      if (InterpolatedField[rkph] == NULL) {
	ENZO_VFAIL("InterpolatedField[%"ISYM"] not allocated.\n", rkph)
      }
      if (InterpolatedField[rgamma] == NULL) {
	ENZO_VFAIL("InterpolatedField[%"ISYM"] not allocated.\n", rgamma)
      }

      if (rkph != FieldUndefined)
	this->ComputeCellCenteredField(rkph);

      if (rgamma != FieldUndefined)
	this->ComputeCellCenteredField(rgamma);

    } // ENDIF fields to interpolate

  /* Clean up :: these are defined in ComputeCellCenteredField */

  if (BaryonField[NumberOfBaryonFields] != NULL) {

    delete [] BaryonField[NumberOfBaryonFields];
    delete [] InterpolatedField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
    InterpolatedField[NumberOfBaryonFields] = NULL;
  }
  
  return SUCCESS;

}
