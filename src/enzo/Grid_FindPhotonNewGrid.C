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

int grid::FindPhotonNewGrid(grid **Grids0, int nGrids0, FLOAT *r, 
			    const FLOAT *u, PhotonPackageEntry* &PP,
			    grid* &MoveToGrid, int &DeltaLevel,
			    const float *DomainWidth, int &DeleteMe,
			    grid *ParentGrid)
{

  int dim, dummy;

  switch (GravityBoundaryType) {
  case TopGridPeriodic: 
    FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], u[0], u[1], u[2]);
    if (dummy >= 0) {
      MoveToGrid = Grids0[dummy];
      DeltaLevel = 0;
      PP->Radius += PFLOAT_EPSILON;
      // Sometimes a photon can loiter on the grid boundary until kicked off
      if (MoveToGrid != this)
	return FALSE;
      else
	MoveToGrid = NULL;
    } else
      // Wrap the photon around the boundary
      if (RadiativeTransferPeriodicBoundary) {
	for (dim = 0; dim < 3; dim++) {
	  if (r[dim] < DomainLeftEdge[dim]) {
	    PP->SourcePosition[dim] += DomainWidth[dim];
	    r[dim] += DomainWidth[dim];
	  } else if (r[dim] > DomainRightEdge[dim]) {
	    PP->SourcePosition[dim] -= DomainWidth[dim];
	    r[dim] -= DomainWidth[dim];
	  }
	} // ENDFOR dim
	FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], 
		     u[0], u[1], u[2]);
	if (dummy == INT_UNDEFINED)
	  ENZO_FAIL("Periodic photon boundary failed to find the next grid.");
	MoveToGrid = Grids0[dummy];
	DeltaLevel = 0;
	PP->Radius += PFLOAT_EPSILON;
	// Sometimes a photon can loiter on the grid boundary until kicked off
	if (MoveToGrid != this)
	  return FALSE;
	else
	  MoveToGrid = NULL;
      } // ENDIF periodic boundary

      else {
	// PhotonPackage left the box
	PP->Photons=-1;
	DeleteMe = TRUE;
	return FALSE;
      }

  case TopGridIsolated: 
    FindRootGrid(dummy, Grids0, nGrids0, r[0], r[1], r[2], 
		 u[0], u[1], u[2]);
    if (dummy >= 0) {
      MoveToGrid = Grids0[dummy];
      DeltaLevel = 0;
      PP->Radius += PFLOAT_EPSILON;
    } else {
      // PhotonPackage left the box
      PP->Photons=-1;
      DeleteMe = TRUE;
    }
    if (MoveToGrid != this)
      return FALSE;
    else
      MoveToGrid = NULL;
      
  case SubGridIsolated:
    MoveToGrid = ParentGrid;
    DeltaLevel = -1;
    PP->Radius += PFLOAT_EPSILON;
    if (DEBUG) 
      fprintf(stdout, "Walk: left grid: sent photon to grid %x\n", 
	      ParentGrid);
    return FALSE;

  case GravityUndefined:
  default:
    fprintf(stdout, "grid::WalkPhotonPackage: "
	    "GravityBoundaryType = RadiationBoundary undefined %"ISYM".\n",
	    GravityBoundaryType);
    ENZO_FAIL("");
  } // ENDSWITCH

  return TRUE;

}
