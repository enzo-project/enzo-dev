#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (DELETE EMISSIVITY FIELDS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/  PURPOSE: 
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

int FindField(int f, int farray[], int n);

int grid::DeleteEmissivity(void)
{

  int field, FieldNum, FirstField, LastField;

  if (RadiativeTransfer == FALSE || RadiativeTransferFLD == FALSE) 
    return SUCCESS;

  FirstField = Emissivity0;
  LastField = Emissivity0;

  for (field = FirstField; field < LastField; field++) {

    FieldNum = FindField(field, FieldType, NumberOfBaryonFields);
    if (MyProcessorNumber == ProcessorNumber) {
      delete [] BaryonField[FieldNum];
      BaryonField[FieldNum] = NULL;
    }

    FieldType[FieldNum] = FieldUndefined;

  } // ENDFOR field

  NumberOfBaryonFields -= LastField - FirstField + 1;
	 
  return SUCCESS;
}
