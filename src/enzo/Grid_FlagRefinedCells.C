/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Rick Wagner                                                *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Board of Trustees of the University of Illinois            *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS (FLAG CELLS THAT OVERLAP REFINED CELLS)
/
/  written by: Rick Wagner
/  date:       November, 2005
/  modified1:
/
/  PURPOSE:
/
/
************************************************************************/

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

int grid::FlagRefinedCells(grid *Subgrid)
{

  /* Return if this grid is not on this processor,
   or if there is no subgrid. */

  if (MyProcessorNumber != ProcessorNumber 
      || Subgrid == NULL )
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index;
  int SubgridStart[MAX_DIMENSION], SubgridEnd[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    SubgridStart[dim] = 0;
    SubgridEnd[dim] = 0;
  }

  /* Compute start and stop indices of the active region of the subgrid
     within this grid (and check to make sure subgrid is within this grid). */

  for (dim = 0; dim < GridRank; dim++) {

    if (Subgrid->GridRightEdge[dim] <= GridLeftEdge[dim] ||
	Subgrid->GridLeftEdge[dim]  >= GridRightEdge[dim])
      return SUCCESS;

    SubgridStart[dim] = nint(
        (Subgrid->GridLeftEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim];
    SubgridEnd[dim] = nint(
	(Subgrid->GridRightEdge[dim] - GridLeftEdge[dim])/CellWidth[dim][0]
			       ) + GridStartIndex[dim] - 1;

    SubgridStart[dim] = max(SubgridStart[dim], GridStartIndex[dim]);
    SubgridEnd[dim]   = min(SubgridEnd[dim], GridEndIndex[dim]);

  }

  /* Now that there is overlap, set the flagging field to true. */
  
  for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
    for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	FlaggingField[index + i] = 1;
    }

  return SUCCESS;
  
}
