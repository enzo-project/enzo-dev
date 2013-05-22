/***********************************************************************
/
/  GRID CLASS (FLAGS CELL FOR REFINEMENT DEPENDING ON ITS REGION)
/
/  written by: Elizabeth Tasker
/  date:       May, 2013
/  modified1:
/
/  PURPOSE: Allows for non-cubic geometries in refined region at different levels
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

  /* Return if this grid is not on this processor
     or if multi refined regions are not being used . */

  if (MyProcessorNumber != ProcessorNumber || MultiRefineRegionGeometry[0] < 0)
    return SUCCESS;


  /* declarations */

  int i, j, k, index, dim, region;
  float xpos, ypos, zpos, ring_width, rad, dr;
  int LocalMaximumRefinementLevel = 0;


  /* Default values */

  if (MultiRefineRegionMaximumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMaximumOuterLevel = MaximumRefinementLevel;
  if (MultiRefineRegionMinimumOuterLevel == INT_UNDEFINED)
    MultiRefineRegionMinimumOuterLevel = 0;


  /* error check */

  if (FlaggingField == NULL) 
    ENZO_FAIL("Flagging Field is undefined");

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Create temporary array for cell flagging on this grid */

  int *TempFlaggingField = new int[size];
  for (i = 0; i < size; i++)
    TempFlaggingField[i] = -1;

  /* loop over multiregions */
  
  for (region = 0; region < MAX_STATIC_REGIONS; region++) {

    if (MultiRefineRegionMaximumLevel[region] == INT_UNDEFINED)
      MultiRefineRegionMaximumLevel[region] = MaximumRefinementLevel;

    /* loop over grid */

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

	index = (k*GridDimension[1] + j)*GridDimension[0] + 
	  GridStartIndex[0];

	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	   switch (MultiRefineRegionGeometry[region]) { 

	   case 0: // cubic geometry. Default value.

	     if (GridRank == 1)
	      if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] > 
		  MultiRefineRegionLeftEdge[region][0]                 &&
		  CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < 
		  MultiRefineRegionRightEdge[region][0]) {
		if (level >= MultiRefineRegionMaximumLevel[region] || level < MultiRefineRegionMinimumLevel[region])
		  if (level < MultiRefineRegionMinimumLevel[region])
		    TempFlaggingField[index] = 1; // flag
		  else
		    TempFlaggingField[index] = 0; // unflag
		else 
		  TempFlaggingField[index] = 2; // leave as it
	      }

	    if (GridRank == 2)
	      if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] > 
		  MultiRefineRegionLeftEdge[region][0]                 &&
		  CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < 
		  MultiRefineRegionRightEdge[region][0]               &&
	      
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] > 
		  MultiRefineRegionLeftEdge[region][1]                 &&
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] < 
		  MultiRefineRegionRightEdge[region][1]) {
		if (level >= MultiRefineRegionMaximumLevel[region] || level < MultiRefineRegionMinimumLevel[region])
		  if (level < MultiRefineRegionMinimumLevel[region])
		    TempFlaggingField[index] = 1; // flag
		  else
		    TempFlaggingField[index] = 0; // unflag
		else 
		  TempFlaggingField[index] = 2; // leave as it
	      }

	    if (GridRank > 2) 
	      if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] > 
		  MultiRefineRegionLeftEdge[region][0]                 &&
		  CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < 
		  MultiRefineRegionRightEdge[region][0]               &&
	      
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] > 
		  MultiRefineRegionLeftEdge[region][1]                 &&
		  CellLeftEdge[1][j] + 0.5*CellWidth[1][j] < 
		  MultiRefineRegionRightEdge[region][1]               && 
		
		  CellLeftEdge[2][k] + 0.5*CellWidth[2][k] > 
		  MultiRefineRegionLeftEdge[region][2]                &&
		  CellLeftEdge[2][k] + 0.5*CellWidth[2][k] < 
		  MultiRefineRegionRightEdge[region][2]){
		if (level >= MultiRefineRegionMaximumLevel[region] || level < MultiRefineRegionMinimumLevel[region])
		  if (level < MultiRefineRegionMinimumLevel[region])
		    TempFlaggingField[index] = 1; // flag
		  else
		    TempFlaggingField[index] = 0; // unflag
		else 
		  TempFlaggingField[index] = 2; // leave as it
	      }
	
 	    break;
	    
	   case 1: // ring
	   case 2: // cylinder 

	     /* Check parameters are set up correctly. Set some defaults */

	     if (GridRank < 3)  // check 3D
	       ENZO_FAIL("MultiRefineRegionGeometry for ring (2) or cyclinder (3)  must be 3D... well OK, it could be 2D but you'll need to adjust the code.");

	     if (MultiRefineRegionOrientation[region][0] == FLOAT_UNDEFINED || MultiRefineRegionRadius[region] == FLOAT_UNDEFINED)
	       ENZO_FAIL("Parameters for refined region not correctly set for geometry. Need MultiRefineRegionOrientation, MultiRefineRegionRadius (MultiRefineRegionWidth for ring.\n");

	     if (MultiRefineRegionCenter[region][0] == FLOAT_UNDEFINED)
	       for (dim = 0; dim < GridRank; dim++)
		 MultiRefineRegionCenter[region][dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
	     
	     /* Calculate position */

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
	     LocalMaximumRefinementLevel = -1;

	     if (MultiRefineRegionGeometry[region] == 1) { 

	       /* Refine gradually inside doughnut from max outside value to inner max so no nasty boundary mess-ups */
	       /* Note, if you don't want the staggering (e.g. for a regions within regions set-up) set 
		  MultiRefineRegionStaggeredRefinement[region] = 0 */

	       dr = fabs(MultiRefineRegionRadius[region]-rad);
	       if (dr < 0.5*MultiRefineRegionWidth[region]) { // inside doughnut

		 if (dr < 0.5*MultiRefineRegionWidth[region]-MultiRefineRegionStaggeredRefinement[region]) // flat region at maximum refinement
		   LocalMaximumRefinementLevel = MultiRefineRegionMaximumLevel[region];
		 else
		   LocalMaximumRefinementLevel = int(MultiRefineRegionMaximumLevel[region] - (dr-(0.5*MultiRefineRegionWidth[region]-MultiRefineRegionStaggeredRefinement[region]))*(MultiRefineRegionMaximumLevel[region]-MultiRefineRegionMaximumOuterLevel)/MultiRefineRegionStaggeredRefinement[region]);
		 
		 if (level >= LocalMaximumRefinementLevel || level < MultiRefineRegionMinimumLevel[region])
		   if (level < MultiRefineRegionMinimumLevel[region])
		     TempFlaggingField[index] = 1; // flag
		   else
		     TempFlaggingField[index] = 0; // unflag
		 else 
		   TempFlaggingField[index] = 2; // leave as it

	       }	     
	     } // end if Geometry = 1

	     else if (MultiRefineRegionGeometry[region] == 2) { // cylinder
	       
	       if (rad < MultiRefineRegionRadius[region]) { // inside cyclinder
	       
		 // middle of region, maximum allowed AMR applies
		 if (rad < MultiRefineRegionRadius[region]-MultiRefineRegionStaggeredRefinement[region])
		   LocalMaximumRefinementLevel =  MultiRefineRegionMaximumLevel[region];
		 else 
		   LocalMaximumRefinementLevel = int(MultiRefineRegionMaximumLevel[region] - (rad-(MultiRefineRegionRadius[region]-MultiRefineRegionStaggeredRefinement[region]))*(MultiRefineRegionMaximumLevel[region]-MultiRefineRegionMaximumOuterLevel)/MultiRefineRegionStaggeredRefinement[region]);

		 if (level >= LocalMaximumRefinementLevel || level < MultiRefineRegionMinimumLevel[region])
		   if (level < MultiRefineRegionMinimumLevel[region])
		     TempFlaggingField[index] = 1; // flag
		   else
		     TempFlaggingField[index] = 0; // unflag
		 else 
		   TempFlaggingField[index] = 2; // leave as it

	       }
	     } // end if Geometry = 2

	     break;

	   default:
	     if (fabs(MultiRefineRegionGeometry[region]) < 0) 
	       ENZO_FAIL("MultiRefineRegionGeometry undefined");
	       
	   } // end switch

	} 
      } // end loop over grid
  } // end loop over regions

  
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {

    if (TempFlaggingField[i] < 0) { // outside all regions
      
      if (level >= MultiRefineRegionMaximumOuterLevel)
	FlaggingField[i] = 0;
      else if (level < MultiRefineRegionMinimumOuterLevel)
	FlaggingField[i]++;
    }
    else if (TempFlaggingField[i] == 0)
      FlaggingField[i] = 0;
    else if (TempFlaggingField[i] == 1)
      FlaggingField[i] = 1;
    else 
      FlaggingField[i] = FlaggingField[i]; // leave as is

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


  delete [] TempFlaggingField;

  return SUCCESS;

}
