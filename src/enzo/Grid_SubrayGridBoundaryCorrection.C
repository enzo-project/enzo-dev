#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CORRECTION TO SUBRAY RADIUS AT GRID BOUNDARIES)
/
/  written by: John Wise
/  date:       September, 2010
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::SubrayGridBoundaryCorrection
(PhotonPackageEntry* &PP, FLOAT u[], int u_sign[], int u_dir[], 
 FLOAT u_inv[], FLOAT f[], FLOAT r[], int g[])
{

  int dim;
  bool any_face = false;
  bool on_left_face[] = {false, false, false};
  bool on_right_face[] = {false, false, false};

  for (dim = 0; dim < GridRank; dim++) {
    if (g[dim] <= GridStartIndex[dim])
      on_left_face[dim] = true;
    else if (g[dim] >= GridEndIndex[dim])
      on_right_face[dim] = true;
    if (on_left_face[dim] || on_right_face[dim])
      any_face = true;
  }

  /* If the ray isn't on a cell boundary but is on a face cell, assume
     it has been recently created from a splitting main ray. */

  bool inside_cell = true;
  for (dim = 0; dim < GridRank; dim++)
    inside_cell &= !(r[dim] == f[dim]);

  if (any_face && inside_cell && PP->Radius > 0) {

    int direction, u_ndir[3];
    FLOAT lce[3], dri[3];
    FLOAT max_dr;

    // Determine last cell boundary crossing (! negates the direction
    // of the ray)
    for (dim = 0; dim < 3; dim++)
      lce[dim] = CellLeftEdge[dim][g[dim] - !u_dir[dim]];

    // Radius of the last edge crossing in each dimension
    for (dim = 0; dim < 3; dim++)
      dri[dim] = u_inv[dim] * (lce[dim] - PP->SourcePosition[dim]);

    // The closest one is the one we want
    for (dim = 1, direction = 0, max_dr = dri[0]; dim < 3; dim++)
      if (dri[dim] > max_dr) {
	direction = dim;
	max_dr = dri[dim];
      }

    /* If the direction of the closest cell wall and grid wall are the
       same, the main ray was most likely split at the grid boundary.
       We adjust the ray's radius to correspond to the grid
       boundary. */

    if ((on_left_face [direction] && u_dir[direction] == 1) ||
	(on_right_face[direction] && u_dir[direction] == 0))
      PP->Radius = max_dr;

  } // ENDIF inside a cell on the grid face

  return SUCCESS;

}
