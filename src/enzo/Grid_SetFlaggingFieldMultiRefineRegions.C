/***********************************************************************
/
/  GRID CLASS (SETS FLAG TO ZERO IF CELL IS OUTSIDE ALLOWED REFINE REGION)
/
/  written by: Elizabeth Tasker
/  date:       July, 2010
/  modified1:
/
/  PURPOSE: Allows for non-cubic geometries in refined region
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::SetFlaggingFieldMultiRefineRegions(int level)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */

  int i, j, k, index, dim, region;
  float xpos, ypos, zpos, ring_width, rad, dr;
  int RegionMaximumRefinementLevel = MaximumRefinementLevel;
  int NewMaximumRefinementLevel = MaximumRefinementLevel;
  int RegionMinimumRefinementLevel = 0;

  /* error check */

  if (FlaggingField == NULL) 
    ENZO_FAIL("Flagging Field is undefined");
    
  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

 
  /* loop over multiregions */
  
  for (region = 0; region < MAX_STATIC_REGIONS; region++) {

    if (MultiRefineRegionGeometry[region] < 0)
      continue;

    /* Find maximum level */
    /* The default is the root grid, 0, but if there are static regions, their highest level is used instead */

    if (MultiRefineRegionMaximumLevel[region] > 0) {
      RegionMaximumRefinementLevel = MultiRefineRegionMaximumLevel[region]; 
      NewMaximumRefinementLevel = MultiRefineRegionMaximumLevel[region];
    }
    
    if (MultiRefineRegionMinimumLevel[region] < 0) {
      for (i = 0; i < MAX_STATIC_REGIONS; i++)
	if (StaticRefineRegionLevel[i] > RegionMinimumRefinementLevel)
	  RegionMinimumRefinementLevel = StaticRefineRegionLevel[i];
      RegionMinimumRefinementLevel++; // should be one above outer layer of ring
    }
    else 
      RegionMinimumRefinementLevel = MultiRefineRegionMinimumLevel[region];

    /* loop over grid */

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

	index = (k*GridDimension[1] + j)*GridDimension[0] + 
	  GridStartIndex[0];

	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	  /* Normal cubic geometry. This is the default value. */
	  switch (MultiRefineRegionGeometry[region]) { 

	  case 0: // rectangular box
	  case 1:
	  
	    if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < 
		MultiRefineRegionLeftEdge[region][0]                  ||
		CellLeftEdge[0][i] + 0.5*CellWidth[0][i] > 
		MultiRefineRegionRightEdge[region][0]) 
	      FlaggingField[index] = 0;

	    if (GridRank > 1)
	      if (CellLeftEdge[1][j] + 0.5*CellWidth[1][j] < 
		  MultiRefineRegionLeftEdge[region][1]                  ||
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] > 
		  MultiRefineRegionRightEdge[region][1]) 
		FlaggingField[index] = 0;

	    if (GridRank > 2)
	      if (CellLeftEdge[2][k] + 0.5*CellWidth[2][k] < 
		  MultiRefineRegionLeftEdge[region][2]                  ||
		  CellLeftEdge[2][k] + 0.5*CellWidth[2][k] > 
		  MultiRefineRegionRightEdge[region][2]) 
		FlaggingField[index] = 0;
	  
	    break;

	  case -1: // avoid refinement box

	    if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] > 
		MultiRefineRegionLeftEdge[region][0]                  ||
		CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < 
		MultiRefineRegionRightEdge[region][0]) 
	      FlaggingField[index] = 0;

	    if (GridRank > 1)
	      if (CellLeftEdge[1][j] + 0.5*CellWidth[1][j] > 
		  MultiRefineRegionLeftEdge[region][1]                  ||
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] < 
		  MultiRefineRegionRightEdge[region][1]) 
		FlaggingField[index] = 0;

	    if (GridRank > 2)
	      if (CellLeftEdge[2][k] + 0.5*CellWidth[2][k] > 
		  MultiRefineRegionLeftEdge[region][2]                  ||
		  CellLeftEdge[2][k] + 0.5*CellWidth[2][k] < 
		  MultiRefineRegionRightEdge[region][2]) 
		FlaggingField[index] = 0;
	  
	    break;

       	
	    /* 2, 3 cases: cyclindrical geometry for which disk radius is needed */
	  case 2: // ring 
	  case -3: // cylinder in which there is NO refinement

	    if (GridRank < 3)  // check 3D
	      ENZO_FAIL("MultiRefineRegionGeometry for ring (2) must be 3D... well OK, it could be 2D but you'll need to adjust the code."); 

	    if (MultiRefineRegionOrientation[region][0] == FLOAT_UNDEFINED) 
	      ENZO_FAIL("Need MultiRefineRegionOrientation for orientation in MultiRefineRegionGeometry = 2+\n");

	    	 
	    xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - 
	      MultiRefineRegionCenter[region][0];
	  
	    ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -
	      MultiRefineRegionCenter[region][1];
		
	    zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - 
	      MultiRefineRegionCenter[region][2];

	 
	    /* Compute distance from center in disk plane */
	    xpos -= xpos*MultiRefineRegionOrientation[region][0];
	    ypos -= ypos*MultiRefineRegionOrientation[region][1];
	    zpos -= zpos*MultiRefineRegionOrientation[region][2];

	    rad = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
	  	  	  
	    
	    if (MultiRefineRegionGeometry[region] == 2) { 

	      if ((rad < MultiRefineRegionRadius[region] - 0.5*MultiRefineRegionWidth[region]) || 
		  (rad > MultiRefineRegionRadius[region] + 0.5*MultiRefineRegionWidth[region]) ||
		  fabs(zpos) > 1.0)
		FlaggingField[index] = 0;
	  
	      /* Refine gradually inside doughnut so no nasty boundary mess-ups */
	  
	      dr = fabs(MultiRefineRegionRadius[region]-rad);
	    if (dr < 0.5*MultiRefineRegionWidth[region]-MultiRefineRegionStaggeredRefinement[region]) // flat region at maximum refinement
	      RegionMaximumRefinementLevel = NewMaximumRefinementLevel;
	    else
	      RegionMaximumRefinementLevel = int(NewMaximumRefinementLevel - (dr-(0.5*MultiRefineRegionWidth[region]-MultiRefineRegionStaggeredRefinement[region]))*(NewMaximumRefinementLevel-RegionMinimumRefinementLevel)/MultiRefineRegionStaggeredRefinement[region]);
		
	   
	    } // end if (MultiRefineRegionGeometry == 2)
	    else if (MultiRefineRegionGeometry[region] == -3) {

	      if (rad < MultiRefineRegionRadius[region]-MultiRefineRegionStaggeredRefinement[region])
		FlaggingField[index] = 0;

	      /* Refine gradually upto inner no-refine doughnut so no nasty boundary mess-ups */
	      if (rad > MultiRefineRegionRadius[region])
		RegionMaximumRefinementLevel = NewMaximumRefinementLevel;
	      else
		RegionMaximumRefinementLevel = int(RegionMinimumRefinementLevel + (rad-(MultiRefineRegionRadius[region]-MultiRefineRegionStaggeredRefinement[region]))*(NewMaximumRefinementLevel-RegionMinimumRefinementLevel)/MultiRefineRegionStaggeredRefinement[region]);
	    } // end if (MultiRefineRegionGeometry == -3)

	  
	    if (level >= RegionMaximumRefinementLevel)
	      FlaggingField[index] = 0;
	  
	    break;

	    
	  default:
	    ENZO_FAIL("MultiRefineRegionGeometry undefined");

	  } // end switch

	  FlaggingField[index] = min(FlaggingField[index], 1); 

	} 
      } // end loop over grid
  } // end loop over regions
  
  /* If debuging, count up the number of flagged cells & report. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    if (FlaggingField[i] > 1)
      FlaggingField[i] = 1;
    NumberOfFlaggedCells += FlaggingField[i];
  }
  
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridEndIndex[dim] - GridStartIndex[dim] + 1;

  if (debug)
    printf("SetFlaggingFieldMultiRefineRegions: NumberOfFlaggedCells = %d (%.1f%%)\n",
	   NumberOfFlaggedCells, float(NumberOfFlaggedCells)*100.0/
	   float(size));


  return SUCCESS;

}
