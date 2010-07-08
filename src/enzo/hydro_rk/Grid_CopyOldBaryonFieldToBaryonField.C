#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::CopyOldBaryonFieldToBaryonField()
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int size = 1;
  for (int dim = 0; dim < GridDimension[dim]; dim++)
    size *= GridDimension[dim];

  for (int field = 0; field < NumberOfBaryonFields; field++) {

    if (OldBaryonField[field] == NULL) {
      fprintf(stderr, "OldBaryonField missing.\n");
      return FAIL;
    }

    for (int i = 0; i < size; i++)
      BaryonField[field][i] = OldBaryonField[field][i];

  }

  return SUCCESS;

}

