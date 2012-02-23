/***********************************************************************
/
/  GRID: ADD RADIO-MODE JET-LIKE FEEDBACK BASED ON STATIC SMBH
/
/  written by: Yuan Li and Greg Bryan
/  date:       December, 2011
/  modified1: 
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"


int grid::ClusterSMBHFeedback(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Return if not on most-refined level. */

  if (level != MaximumRefnementLevel)
    return SUCCESS;

  /* Compute the jet launching region
     (assume jet launched from PointSourceGravityPosition) */

  FLOAT JetLeftCorner[MAX_DIMENSION], LeftRightCorner[MAX_DIMENSION];
  FLOAT JetCenter[MAX_DIMENSION];
  int jet_dim = 2;  // z-axis (should make parameter?)
  int JetLaunchOffset = 10; // 10 cellwidth
  int JetRadius = 3; // cellwidths

  for (dim = 0; dim < GridRank; dim++) {
    JetCenter[dim] = PointSourceGravityPosition[dim];
    JetLeftCorner[dim] = JetCenter[dim];
    JetRightCorner[dim] = JetCenter[dim];
    if (dim != jet_dim) {
      JetLeftCorner[dim] -= JetRadius*CellWidth[dim][0];
      JetRightCorner[dim] += JetRadius*CellWidth[dim][0];
    }
  }

  JetLeftCorner[jet_dim] -= JetLaunchOffset*CellWidth[jet_dim][0];
  JetRightCorner[jet_dim] += JetLaunchOffset*CellWidth[jet_dim][0];
 
  /* Compute indices of jet launch region. */

  int JetStartIndex[MAX_DIMENSION], JetEndIndex[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {

    /* Compute start and end indices of jet */

    JetStartIndex[dim] = nint((JetLeftCorner[dim] - CellLeftEdge[dim][0])/CellWidth[dim][0]);
    JetEndIndex[dim] = nint((JetRightCorner[dim] - CellLeftEdge[dim][0])/CellWidth[dim][0]);

    /* If Jet is not on this grid, return. */

    if (JetStartIndex[dim] > GridDimension[dim]-1 || JetEndIndex[dim] < 0)
      return SUCCESS;

    /* Clip edge of jet launching disk */

    if (dim != jet_dim) {
      JetStartIndex[dim] = max(JetStartIndex[dim], 0);
      JetEndIndex[dim] = min(JetStartIndex[dim], GridDimension[dim]-1);
    }

  } // end: loop over dim

  /* Compute mass and momentum to be put into cells in code units. */

  

  /* Loop over launch disks and set cell values (this code assumes jet_dim = 2). */

  if (JetStartIndex[jet_dim] >= 0) {
    k = JetStartIndex[jet_dim];
    for (j = JetStartIndex[1]; j <= JetEndIndex[1]; j++) {
      for (i = JetStartIndex[0]; i <= JetEndIndex[0]; i++) {
	radius = sqrt(pow((CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - JetCenter[0]), 2) + 
		      pow((CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - JetCenter[1]), 2))
	  /CellWidth[0][0]; // in cell widths
	
	BaryonField[DensNum][GRIDINDEX(i,j,k)] += XXX*exp(radius/JetRadius);
	BaryonField[Vel3Num][GRIDINDEX(i,j,k)] += XXX;
	//	BaryonField[GENum][GRIDINDEX(i,j,k)] += XXX;
      }
    }
  }

    

  /* loop over cells to be modified, add jet mass, momentum, and energy. */

  return SUCCESS;

}
 
