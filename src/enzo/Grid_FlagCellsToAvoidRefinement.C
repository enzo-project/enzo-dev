/***********************************************************************
/
/  GRID CLASS (FLAG CELLS IN WHICH TO AVOID REFINEMENT)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int FindField(int field, int farray[], int numfields);

int grid::FlagCellsToAvoidRefinement()
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */

  if (FlaggingField == NULL) {
    ENZO_FAIL("Flagging field is undefined.");
  }

  /* compute size */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  int NumberOfFlaggedCells = 0;
  int ColorField = FindField(ForbiddenRefinement, FieldType, NumberOfBaryonFields); 
  if (ColorField == -1)
    ENZO_FAIL("Can't avoid refinement without ForbiddenRefinement field!");
  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields); 
  if (ColorField == -1)
    ENZO_FAIL("Can't find Density field.");
 
  for (i = 0; i < size; i++)
    if(BaryonField[ColorField][i] > 0.125 * BaryonField[DensNum][i]) {
      FlaggingField[i] = 0.0;
      //NumberOfFlaggedCells++;
    }
 
  return NumberOfFlaggedCells;
 
}
