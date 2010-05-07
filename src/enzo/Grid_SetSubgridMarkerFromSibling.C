/***********************************************************************
/
/  GRID CLASS (MARK SUBGRID)
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:  John Wise, May, 2010 -- modified the subgrid version to
/              accept sibling grids with a possible periodic offset.
/
/  PURPOSE:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int grid::SetSubgridMarkerFromSibling(grid *Sibling, 
				      FLOAT EdgeOffset[MAX_DIMENSION])
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index, size;
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];

  /* if the field has not been allocated yet, do it here */ 

  if (SubgridMarker == NULL)  {
    for (dim = 0, size = 1; dim < GridRank; dim++)
      size *= GridDimension[dim];
    SubgridMarker = new grid*[size];
    for (i=0; i<size; i++) SubgridMarker[i] = NULL;
  }

  /* Compute the left and right edges of this grid (including ghost zones). */
 
  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = CellLeftEdge[dim][0] + EdgeOffset[dim];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
      CellWidth[dim][GridDimension[dim]-1] +
      EdgeOffset[dim];
  }

  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= Sibling->GridRightEdge[dim] ||
        GridRight[dim] <= Sibling->GridLeftEdge[dim]   )   
      return SUCCESS;

  /* Compute start and stop indices of the active region of the subgrid
     within this grid (and check to make sure subgrid is within this grid). */

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute left and right positions in problem space.  note:
       include buffer zones of this grid but not the other grid. */
 
    Left[dim]  = max(GridLeft[dim], Sibling->GridLeftEdge[dim]);
    Right[dim] = min(GridRight[dim], Sibling->GridRightEdge[dim]);

    /* Convert this to index positions in this grid */
 
    Start[dim] = nint((Left[dim]  - GridLeft[dim]) / CellWidth[dim][0]);
    End[dim]   = nint((Right[dim] - GridLeft[dim]) / CellWidth[dim][0]) - 1;

    if (End[dim] - Start[dim] < 0)
      return SUCCESS;

  }

  /* Now that there is overlap mark this grid with the subgrid pointer. */

  for (k = Start[2]; k <= End[2]; k++)
    for (j = Start[1]; j <= End[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = Start[0]; i <= End[0]; i++)
	SubgridMarker[index + i] = Sibling;
    }
  
  return SUCCESS;
  
}
