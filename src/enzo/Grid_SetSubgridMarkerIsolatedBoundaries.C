/***********************************************************************
/
/  GRID CLASS (MARK SUBGRID FOR ISOLATED BOUNDARIES)
/
/  written by: John Wise
/  date:       April, 2012
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

int grid::SetSubgridMarkerIsolatedBoundaries(void)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  int i, j, k, n, dim, field, index;
  int SubgridStart[MAX_DIMENSION], SubgridEnd[MAX_DIMENSION];

  /* Check if the outer cells lie outside the domain */

  bool inside = true;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    inside &= CellLeftEdge[dim][0] > DomainLeftEdge[dim] &&
      CellLeftEdge[dim][GridDimension[dim]-1] < DomainRightEdge[dim];
  if (inside) return SUCCESS;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    SubgridStart[dim] = 0;
    SubgridEnd[dim] = 0;
  }

  /* Compute start and stop indices of the domain within this grid
     (and check to make sure subgrid is within this grid). */

  for (dim = 0; dim < GridRank; dim++) {

    SubgridStart[dim] = nint(
        (DomainLeftEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]) 
      + GridStartIndex[dim];
    SubgridEnd[dim] = nint(
	(DomainRightEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim] - 1;

    SubgridStart[dim] = max(SubgridStart[dim], 0);
    SubgridEnd[dim]   = min(SubgridEnd[dim], GridDimension[dim]-1);

  }

  /* Now that there is overlap mark this grid with a NULL pointer. */

  int xdim, ydim, din;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    xdim = (dim+1) % MAX_DIMENSION;
    ydim = (dim+2) % MAX_DIMENSION;

    // minus face
    for (k = 0; k < SubgridStart[dim]; k++)
      for (j = 0; j < GridDimension[ydim]; j++)
	for (i = 0; i < GridDimension[xdim]; i++) {
	  switch (dim) {
	  case 0:
	    index = GRIDINDEX_NOGHOST(k,i,j);
	    break;
	  case 1:
	    index = GRIDINDEX_NOGHOST(j,k,i);
	    break;
	  case 2:
	    index = GRIDINDEX_NOGHOST(i,j,k);
	    break;
	  }
	  SubgridMarker[index] = NULL;
	} // ENDFOR i
    
    // plus face
    for (k = SubgridEnd[dim]+1; k < GridDimension[dim]; k++)
      for (j = 0; j < GridDimension[ydim]; j++)
	for (i = 0; i < GridDimension[xdim]; i++) {
	  switch (dim) {
	  case 0:
	    index = GRIDINDEX_NOGHOST(k,i,j);
	    break;
	  case 1:
	    index = GRIDINDEX_NOGHOST(j,k,i);
	    break;
	  case 2:
	    index = GRIDINDEX_NOGHOST(i,j,k);
	    break;
	  }
	  SubgridMarker[index] = NULL;
	} // ENDFOR i

  } // ENDFOR dim

  return SUCCESS;
  
}
