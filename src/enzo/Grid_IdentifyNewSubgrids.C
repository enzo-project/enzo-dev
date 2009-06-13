/***********************************************************************
/
/  GRID CLASS (IDENTIFY A LIST OF POSSIBLE SUBGRIDS, SIMPLE VERSION)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
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
 
void grid::IdentifyNewSubgrids(GridList &SubGridList)
{
  /* declarations */
 
  int i, j, k, dim, InsideRegion, Offset;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  GridList *NewSubgrid;
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return;
 
  /* loop over only one zone inside boundary */
 
  for (dim = 0; dim < 3; dim++) {
    Start[dim] = GridStartIndex[dim];
    End[dim] = GridEndIndex[dim];
  }
 
  /* loop over grid and for every flagged cell, create a new Subgrid */
 
    for (k = Start[2]; k <= End[2]; k++)
      for (j = Start[1]; j <= End[1]; j++)
	for (i = Start[0]; i <= End[0]; i++)
	  if (*(FlaggingField + i + j*GridDimension[0] +
		k*GridDimension[1]*GridDimension[0])      ) {
 
	    /* check if we are inside the allowed refining region. */
 
	    InsideRegion = TRUE;
	    for (dim = 0; dim < GridRank; dim++)
	      if (GridDimension[dim] > 1) {
	       Offset = (dim == 0) ? i : j;
	       Offset = (dim == 2) ? k : Offset;
	       if (*(CellLeftEdge[dim] + Offset) < RefineRegionLeftEdge[dim])
		 InsideRegion = FALSE;
	       if (*(CellLeftEdge[dim] + Offset) >= RefineRegionRightEdge[dim])
		 InsideRegion = FALSE;
	      }
 
	    /* create new subgrid unless we are outside of the region. */
 
	    if (InsideRegion == TRUE) {
 
	      NewSubgrid = new GridList;
 
	      /* set its rank and dimensions */
 
	      NewSubgrid->GridRank = GridRank;
	      for (dim = 0; dim < GridRank; dim++)
		NewSubgrid->GridDimension[dim] = min(GridDimension[dim], 3);
	      NewSubgrid->StartIndex[0] = max(i-1, Start[0]);
	      NewSubgrid->EndIndex[0]   = min(i+1, End[0]);
	      NewSubgrid->StartIndex[1] = max(j-1, Start[1]);
	      NewSubgrid->EndIndex[1]   = min(j+1, End[1]);
	      NewSubgrid->StartIndex[2] = max(k-1, Start[2]);
	      NewSubgrid->EndIndex[2]   = min(k+1, End[2]);
	      NewSubgrid->NumberFlagged = 1;
 
	      /* insert NewSubgrid at head of list */
	
	      NewSubgrid->NextGrid = SubGridList.NextGrid;
	      SubGridList.NextGrid = NewSubgrid;
 
	    }
	    else
	      if (debug)
		fprintf(stderr, "Subgrid rejected due to OutsideRegion.\n");
	
	  }
 
}
