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

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    ENZO_FAIL("Error in: "__FILE__);
  }

  int rkph, rgamma, field;

  for (field = 0; field < NumberOfBaryonFields; field++)
    if (FieldsToInterpolate[field] == TRUE) {
      switch (FieldType[field]) {
      case HIDensity:
	rkph = kphHINum;
	rgamma = gammaHINum;
	break;
      case HeIDensity:
	rkph = kphHeINum;
	rgamma = gammaHeINum;
	break;
      case HeIIDensity:
	rkph = kphHeIINum;
	rgamma = gammaHeIINum;
	break;
      } // ENDSWITCH field

      if (InterpolatedField[rkph] == NULL) {
	fprintf(stderr, "InterpolatedField[%"ISYM"] not allocated.\n", rkph);
	ENZO_FAIL("Error in: "__FILE__);
      }
      if (InterpolatedField[rgamma] == NULL) {
	fprintf(stderr, "InterpolatedField[%"ISYM"] not allocated.\n", rgamma);
	ENZO_FAIL("Error in: "__FILE__);
      }

      if (this->ComputeCellCenteredField(rkph) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeCellCenteredField.\n");
	ENZO_FAIL("Error in: "__FILE__);
      }

      if (this->ComputeCellCenteredField(rgamma) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeCellCenteredField.\n");
	ENZO_FAIL("Error in: "__FILE__);
      }

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
