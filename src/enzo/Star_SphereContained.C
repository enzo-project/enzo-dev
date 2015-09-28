/***********************************************************************
/
/  FIND ACCRETION SPHERE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: 
/
/  PURPOSE: When we remove baryons from the grid to add to the star
/           particle, look for a sphere that contains twice its mass.
/           Stepping outward by a cell width.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"

int Star::SphereContained(LevelHierarchyEntry *LevelArray[], int level, 
			  float Radius)
{

  LevelHierarchyEntry *Temp;
  int i, dim, direction, cornersContained, Rank, result;
  bool inside;
  int cornerDone[8], Dims[MAX_DIMENSION];
  FLOAT corners[MAX_DIMENSION][8];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  /**************************************************************

     Compute corners of cube that contains a sphere of r=Radius

   **************************************************************/

  LCAPERF_START("star_SphereContained");

  for (i = 0; i < 8; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {

      // If the bit is true, forward.  If not, reverse.
      direction = (i >> dim & 1) ? 1 : -1;
      corners[dim][i] = pos[dim] + direction * Radius;
    }
    cornerDone[i] = 0;
  }

  /* Check if the influenced sphere is contained within the grids on
     this level */

  for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {

      for (i = 0; i < 8; i++) {
	if (cornerDone[i]) continue;  // Skip if already locally found
	inside = true;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  inside &= (corners[dim][i] >= Temp->GridData->GridLeftEdge[dim] &&
		     corners[dim][i] <= Temp->GridData->GridRightEdge[dim]);
	if (inside)
	  cornerDone[i] = 1;
      } // ENDFOR corners

  } // ENDFOR grids

  /* Take the MPI_MAX of cornerDone flags, then sum them to see if
     they equal 8.  If so, the sphere is contained within grids on
     this level. */

#ifdef USE_MPI
  //  CommunicationAllReduceValues(cornerDone, 8, MPI_MAX);
#endif

  cornersContained = 0;
  for (i = 0; i < 8; i++)
    cornersContained += cornerDone[i];

  result = (cornersContained == 8);
  LCAPERF_STOP("star_SphereContained");
  return result;

}
