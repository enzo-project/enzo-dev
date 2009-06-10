/***********************************************************************
/
/  GRID CLASS (CLEAR THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// Allocate and clear the flagging field.
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::ClearFlaggingField()
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return;
 
  /* error check */
 
  if (FlaggingField != NULL) {
    fprintf(stderr, "ClearFlaggingField: Warning, field already present.\n");
    delete [] FlaggingField;
  }
 
  /* compute size and allocate */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  FlaggingField = new int[size];
 
  /* Clear it */
 
  for (int i = 0; i < size; i++)
    FlaggingField[i] = 0;
 
}
