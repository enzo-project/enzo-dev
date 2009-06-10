/***********************************************************************
/
/  GRID CLASS (CHECK FOR OVERLAPPING ZONES FROM GRID IN ARGUMENT TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       April, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine checks for overlap between the grid in the argument
//   and the current grid.
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::CopyZonesFromGridCountOnly(grid *OtherGrid, int &Overlap)
{
  /* declarations */
 
  int dim;
  Overlap = FALSE;
 
  /* Compute the left and right edges of this grid (including ghost zones). */
 
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = CellLeftEdge[dim][0];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
                     CellWidth[dim][GridDimension[dim]-1];
  }
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= OtherGrid->GridRightEdge[dim] ||
        GridRight[dim] <= OtherGrid->GridLeftEdge[dim]   )
      return SUCCESS;
 
  /* There is some overlap, so set flag and exit. */
 
  Overlap = TRUE;
 
  return SUCCESS;
}
