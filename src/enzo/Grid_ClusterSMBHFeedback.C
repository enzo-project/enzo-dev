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
  int jet_dim = 2;  // z-axis (should make parameter?)

  for (dim = 0; dim < GridRank; dim++) {
    JetLeftCorner[dim] = PointSourceGravityPosition[dim];
    JetRightCorner[dim] = PointSourceGravityPosition[dim];
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

  /* Loop over launch disks and set cell values. */

  if (JetStartIndex[jet_dim] > 0 && JetStartIndex[jet_dim] < GridDimension[jet_dim]-1) 
    

  /* loop over cells to be modified, add jet mass, momentum, and energy. */

  return SUCCESS;

}
 
