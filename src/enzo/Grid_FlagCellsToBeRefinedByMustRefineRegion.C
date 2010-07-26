/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINE BY REGION)
/
/  written by: Brian O'Shea
/  date:       August 2006
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
 
 
int grid::FlagCellsToBeRefinedByMustRefineRegion(int level)
{
  /* declarations */
 
  int i, j, k, index, dim, size = 1, NumberOfFlaggedCells = 0;
 
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];

  FLOAT CellSize, xpos,ypos,zpos;

  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* error check */
   if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* loop over dimensions - I guess this is unnecessary, 
     but it's handy to have shorter names */
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim] = GridStartIndex[dim];
    End[dim]   = GridEndIndex[dim];
  }


  /* compute size */
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  /* NOTE:  This _only_ works for GridRank = 3 (i.e. 3D calculations).  It could be
     easily extended to 1 and 2D if necessary.  */
  if(GridRank < 3){
    fprintf(stderr, "FlagCellsToBeRefinedByMustRefineRegion only works in 3D!\n");
    return -1;
  }
  
  /* if we're already at the max level that we want to refine to,
     return!  Otherwise we want to refine here. */
  if(MustRefineRegionMinRefinementLevel <= level) return SUCCESS;

  /* check to see if this grid overlaps with our MustRefineRegion. */
  for (dim = 0; dim < GridRank; dim++){

    if(   !((GridRightEdge[dim] > MustRefineRegionLeftEdge[dim]) &&
	    (GridLeftEdge[dim] < MustRefineRegionRightEdge[dim])) )
      return SUCCESS;
    
  }

  /* now we know this grid is at least partially within the refinement region.  
     So we go though each cell, calculate its center, and if it is within the
     refinement region, we flag it! */
  
  CellSize = FLOAT(CellWidth[0][0]);

  /* compute slope */
  for (k = Start[2]; k <= End[2]; k++)
    for (j = Start[1]; j <= End[1]; j++)
      for (i = Start[0]; i <= End[0]; i++) {

	index = i + j*GridDimension[0] +
	  k*GridDimension[1]*GridDimension[0];

	xpos = GridLeftEdge[0] + (FLOAT(i-Start[0])+0.5 )*CellSize;
	ypos = GridLeftEdge[1] + (FLOAT(j-Start[1])+0.5 )*CellSize;
	zpos = GridLeftEdge[2] + (FLOAT(k-Start[2])+0.5 )*CellSize;

	if( (MustRefineRegionLeftEdge[0] <= xpos) && (xpos <= MustRefineRegionRightEdge[0]) &&
	    (MustRefineRegionLeftEdge[1] <= ypos) && (ypos <= MustRefineRegionRightEdge[1]) &&
	    (MustRefineRegionLeftEdge[2] <= zpos) && (zpos <= MustRefineRegionRightEdge[2]) ){

	  FlaggingField[index] += 1;
	}
      }
 
  /* Count number of flagged Cells. */ 
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;
 
}
