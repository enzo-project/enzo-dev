#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CHECKS IF PHOTON NEEDS WRAPPING AROUND PERIODIC BOUNDARY)
/
/  written by: John Wise
/  date:       May, 2010
/  modified1:
/
/  PURPOSE: 
/
/  RETURNS: TRUE  = photon belongs here
/           FALSE = photon should be moved
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

int grid::PhotonPeriodicBoundary(int &cindex, FLOAT *r, int *g, FLOAT *s,
				 PhotonPackageEntry* &PP, grid* &MoveToGrid,
				 const float *DomainWidth, int &DeleteMe)
{

  int dim;
  bool InsideDomain;

  /* Only check for root-level grids. Done outside the routine now. */

//  if (GravityBoundaryType == SubGridIsolated)
//    return TRUE;

  /* Check if the photon is inside the domain */

  for (dim = 0, InsideDomain = true; dim < MAX_DIMENSION; dim++)
    InsideDomain &= (r[dim] >= DomainLeftEdge[dim] && 
		     r[dim] <= DomainRightEdge[dim]);

  /* Return if inside the domain */

  if (InsideDomain)
    return TRUE;

  /* If no periodic boundaries, delete */
  
  if (!RadiativeTransferPeriodicBoundary) {
    PP->Photons = -1;
    MoveToGrid = NULL;
    DeleteMe = TRUE;
    return FALSE;
  }

  /* Now all that's left are photons that have left the domain and
     need to be wrapped around the boundary.  Note that all of these
     photons will remain in the same grid.  This only happens when a
     grid fills the whole domain in one dimension. */

  for (dim = 0; dim < 3; dim++)
    if (r[dim] < DomainLeftEdge[dim]) {
      PP->SourcePosition[dim] += DomainWidth[dim];
      r[dim] += DomainWidth[dim];
      s[dim] += DomainWidth[dim];
      g[dim] += GridEndIndex[dim] - GridStartIndex[dim] + 1;
    } else if (r[dim] > DomainRightEdge[dim]) {
      PP->SourcePosition[dim] -= DomainWidth[dim];
      r[dim] -= DomainWidth[dim];
      s[dim] -= DomainWidth[dim];
      g[dim] -= GridEndIndex[dim] - GridStartIndex[dim] + 1;
    }

  cindex = GRIDINDEX_NOGHOST(g[0], g[1], g[2]);
  MoveToGrid = SubgridMarker[cindex];

  return TRUE;

}
