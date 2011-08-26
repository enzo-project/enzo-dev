/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE AVOID REFINEMENT BY REGION)
/
/  written by: John Wise
/  date:       April, 2011
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
 
 
int grid::FlagCellsToAvoidRefinementRegion(int level)
{
  /* declarations */

  bool outside;
  int i, j, k, index, dim, size, region, NumberOfFlaggedCells = 0;
  int NumberOfRegions = 0;
  FLOAT xpos, ypos, zpos;

  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  for (region = 0; region < MAX_STATIC_REGIONS; region++)
    if (AvoidRefineRegionLevel[region] == INT_UNDEFINED)
      break;
    else
      NumberOfRegions++;

  outside = false;
  for (region = 0; region < NumberOfRegions; region++) {

    if (level >= AvoidRefineRegionLevel[region])
      continue;

    /* check to see if this grid overlaps with this AvoidRefineRegion. */
    for (dim = 0; dim < GridRank; dim++)
      outside |= 
	(!((GridRightEdge[dim] > AvoidRefineRegionLeftEdge[region][dim]) &&
	   (GridLeftEdge[dim] < AvoidRefineRegionRightEdge[region][dim])) );

    if (!outside) break;

  } // ENDFOR regions

  if (outside)
    return SUCCESS;

  /* now we know this grid is at least partially within the refinement region.  
     So we go though each cell, calculate its center, and if it is within the
     avoid refinement region, we unflag it. */

  float half_dx = 0.5 * CellWidth[0][0];
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    zpos = CellLeftEdge[2][k] + half_dx;
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      ypos = CellLeftEdge[1][j] + half_dx;
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	xpos = CellLeftEdge[0][i] + half_dx;

	for (region = 0; region < NumberOfRegions; region++)
	  if( (level >= AvoidRefineRegionLevel[region]) &&
	      (AvoidRefineRegionLeftEdge[region][0] <= xpos) && 
	      (xpos <= AvoidRefineRegionRightEdge[region][0]) &&
	      (AvoidRefineRegionLeftEdge[region][1] <= ypos) && 
	      (ypos <= AvoidRefineRegionRightEdge[region][1]) &&
	      (AvoidRefineRegionLeftEdge[region][2] <= zpos) && 
	      (zpos <= AvoidRefineRegionRightEdge[region][2]) )
	    
	    FlaggingField[index] = 0;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
 
  /* Count number of flagged Cells. */ 
  size = this->GetGridSize();
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;
 
}
