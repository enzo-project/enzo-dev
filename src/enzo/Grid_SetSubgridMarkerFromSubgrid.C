/***********************************************************************
/
/  GRID CLASS (MARK SUBGRID)
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:
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

int grid::SetSubgridMarkerFromSubgrid(grid *Subgrid)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index;
  int SubgridStart[MAX_DIMENSION], SubgridEnd[MAX_DIMENSION];

  /* if the field has not been allocated yet, do it here */ 

  long size = 1;

  for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

  // Subgrid == NULL is the flag to allocate and initialize.
  if (Subgrid == NULL) {
    if (SubgridMarker != NULL)
      delete [] SubgridMarker;
    SubgridMarker = new grid*[size];
    for (i=0; i<size; i++) SubgridMarker[i] = NULL;
    return SUCCESS;
  }

  if (SubgridMarker == NULL)  {
    //    if (debug) printf("allocating SubgridMarker field\n");
    SubgridMarker = new grid*[size];
    if (Subgrid == this) {
      // Set everything to this grid because any cells from subgrids,
      // siblings, or parents will overwrite these.
      for (i=0; i<size; i++) SubgridMarker[i] = this;  //NULL;
//      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
//	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
//	  index = (k*GridDimension[1]+j)*GridDimension[0] + GridStartIndex[0];
//	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
//	    SubgridMarker[index] = this;
//	} // ENDFOR j
      return SUCCESS;

    } else
      for (i=0; i<size; i++) SubgridMarker[i] = NULL;
  }

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    SubgridStart[dim] = 0;
    SubgridEnd[dim] = 0;
  }

  /* Compute start and stop indices of the active region of the subgrid
     within this grid (and check to make sure subgrid is within this grid). */

  for (dim = 0; dim < GridRank; dim++) {

    SubgridStart[dim] = nint(
        (Subgrid->GridLeftEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim];
    SubgridEnd[dim] = nint(
	(Subgrid->GridRightEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim] - 1;

    SubgridStart[dim] = max(SubgridStart[dim], 0);
    SubgridEnd[dim]   = min(SubgridEnd[dim], GridDimension[dim]-1);

  }

  /* Now that there is overlap mark this grid with the subgrid pointer. */

  for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
    for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	SubgridMarker[index + i] = Subgrid;
    }
  
  return SUCCESS;
  
}
