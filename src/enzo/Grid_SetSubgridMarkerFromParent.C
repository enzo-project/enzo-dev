/***********************************************************************
/
/  GRID CLASS (MARK SUBGRID)
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:  John Wise, May, 2010 -- modified the subgrid version to
/              accept parent grids.
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

int grid::SetSubgridMarkerFromParent(grid *Parent)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index, face, side;
  int GZStart[MAX_DIMENSION], GZEnd[MAX_DIMENSION];

  /* if the field has not been allocated yet, do it here */ 

  long size = 1;

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (SubgridMarker == NULL)  {
    SubgridMarker = new grid*[size];
    for (i=0; i<size; i++) SubgridMarker[i] = NULL;
  }

  /* Mark the all ghost zones with the parent grid pointer */

  for (dim = 0; dim < 3; dim++) {

    // Left face
    for (i = 0; i < 3; i++) {
      GZStart[i] = 0;
      GZEnd[i] = (i == dim) ? GridStartIndex[i]-1 : GridDimension[i]-1;
    }

    for (k = GZStart[2]; k <= GZEnd[2]; k++)
      for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = GZStart[0]; i <= GZEnd[0]; i++)
	  SubgridMarker[index + i] = Parent;
      }

    // Right face
    for (i = 0; i < 3; i++) {
      GZStart[i] = (i == dim) ? GridEndIndex[i]+1 : 0;
      GZEnd[i] = GridDimension[i]-1;
    }

    for (k = GZStart[2]; k <= GZEnd[2]; k++)
      for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = GZStart[0]; i <= GZEnd[0]; i++)
	  SubgridMarker[index + i] = Parent;
      }

  } // ENDFOR dim

  /* Now that there is overlap mark this grid with the subgrid pointer. */

  
  return SUCCESS;
  
}
