#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CHECKS IF PHOTON HAS LEFT THE GRID)
/
/  written by: John Wise
/  date:       October, 2009
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
#include "CosmologyParameters.h"
#include "RadiativeTransferHealpixRoutines.h"

int FindRootGrid(int &dummy, grid **Grids0, int nGrids0, FLOAT rx, 
		 FLOAT ry, FLOAT rz, FLOAT ux, FLOAT uy, FLOAT uz);

int grid::FindPhotonNewGrid(int cindex, FLOAT *r, PhotonPackageEntry* &PP,
			    grid* &MoveToGrid, int &DeltaLevel,
			    const float *DomainWidth, int &DeleteMe,
			    grid *ParentGrid)
{

  int dim, dummy, RayInsideGrid;
  bool InsideDomain;

  /* First determine whether the ray has left the grid, store the
     destination grid, and "nudge" the photon package to avoid
     loitering on boundaries. */

  RayInsideGrid = this->PointInGridNB(r);
  MoveToGrid = SubgridMarker[cindex];
  PP->Radius += PFLOAT_EPSILON;

  switch (GravityBoundaryType) {

  case TopGridIsolated: 
  case TopGridPeriodic:
    if (RayInsideGrid) {
      // Inside root grid -> Child grid
      DeltaLevel = +1;
    } else {
      // Outside root grid -> Other root grid or outside domain
      DeltaLevel = 0;
      for (dim = 0, InsideDomain = true; dim < MAX_DIMENSION; dim++)
	InsideDomain &= (r[dim] >= DomainLeftEdge[dim] && 
			 r[dim] <= DomainRightEdge[dim]);
      if (!InsideDomain) {
	if (RadiativeTransferPeriodicBoundary &&
	    GravityBoundaryType == TopGridPeriodic) {
	  for (dim = 0; dim < 3; dim++)
	    if (r[dim] < DomainLeftEdge[dim])
	      PP->SourcePosition[dim] += DomainWidth[dim];
	    else if (r[dim] > DomainRightEdge[dim])
	      PP->SourcePosition[dim] -= DomainWidth[dim];
	} // ENDIF periodic
	else {
	  MoveToGrid = NULL;
	  DeleteMe = TRUE;
	} // ENDELSE periodic
	
      } // ENDIF !InsideDomain
    }
    break;
      
  case SubGridIsolated:
    if (RayInsideGrid) {
      DeltaLevel = +1;
    } else {
      // Outside the grid, we have to determine whether it's a parent
      // or sibling grid
      if (MoveToGrid == ParentGrid)
	DeltaLevel = -1;
      else
	DeltaLevel = 0;
    }
    if (DEBUG) 
      fprintf(stdout, "Walk: left grid: sent photon to grid %x (DeltaL = %d)\n", 
	      MoveToGrid, DeltaLevel);
    break;

  case GravityUndefined:
  default:
    fprintf(stdout, "grid::WalkPhotonPackage: "
	    "GravityBoundaryType = RadiationBoundary undefined %"ISYM".\n",
	    GravityBoundaryType);
    ENZO_FAIL("");
  } // ENDSWITCH

  return SUCCESS;

}
