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

int grid::FindPhotonNewGrid(int cindex, FLOAT *r,
			    PhotonPackageEntry* &PP,
			    grid* &MoveToGrid, int &DeltaLevel,
			    const float *DomainWidth, int &DeleteMe,
			    grid *ParentGrid)
{

  int Refinement;
  int dim, RayInsideGrid;
  bool InsideDomain;

  /* First determine whether the ray has left the grid, store the
     destination grid, and "nudge" the photon package to avoid
     loitering on boundaries. */

  RayInsideGrid = this->PointInGridNB(r);
  MoveToGrid = SubgridMarker[cindex];
  PP->Radius += PFLOAT_EPSILON;

  /* Special case for no gravity (just like TopGridIsolated) */

  if (SelfGravity == FALSE) 
    if (ParentGrid == NULL) {
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
      return TRUE;
    } else {
      // Subgrid
      MoveToGrid = ParentGrid;
      DeltaLevel = -1;
      PP->Radius += PFLOAT_EPSILON;
      if (DEBUG) 
	fprintf(stdout, "Walk: left grid: sent photon to grid %x\n", 
		ParentGrid);
      return FALSE;
    }


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
	if (RadiativeTransferPeriodicBoundary) {
	  for (dim = 0; dim < 3; dim++)
	    if (r[dim] < DomainLeftEdge[dim]) {
 	      PP->SourcePosition[dim] += DomainWidth[dim];
	      //r[dim] += DomainWidth[dim];
	    } else if (r[dim] > DomainRightEdge[dim]) {
 	      PP->SourcePosition[dim] -= DomainWidth[dim];
	      //r[dim] -= DomainWidth[dim];
	    }
	} // ENDIF periodic
	else {
	  MoveToGrid = NULL;
	  DeleteMe = TRUE;
	} // ENDELSE periodic
	
      } // ENDIF !InsideDomain
    } // ENDELSE inside grid
    break;
      
  case SubGridIsolated:
    if (RayInsideGrid) {
      DeltaLevel = +1;
    } else {
      // Outside the grid, we have to determine whether it's a parent
      // or some other grid
      if (MoveToGrid == ParentGrid)
	DeltaLevel = -1;
      else {
	Refinement = nint( MoveToGrid->CellWidth[0][0] /
			   CellWidth[0][0] );
	DeltaLevel = 0;
	while (Refinement > 1) {
	  Refinement /= RefineBy;
	  DeltaLevel--;
	}
      }

      for (dim = 0, InsideDomain = true; dim < MAX_DIMENSION; dim++)
	InsideDomain &= (r[dim] >= DomainLeftEdge[dim] && 
			 r[dim] <= DomainRightEdge[dim]);
      if (!InsideDomain) {
	if (RadiativeTransferPeriodicBoundary) {
	  for (dim = 0; dim < 3; dim++)
	    if (r[dim] < DomainLeftEdge[dim]) {
	      PP->SourcePosition[dim] += DomainWidth[dim];
	      //r[dim] += DomainWidth[dim];
	    } else if (r[dim] > DomainRightEdge[dim]) {
	      PP->SourcePosition[dim] -= DomainWidth[dim];
	      //r[dim] -= DomainWidth[dim];
	    }
	} // ENDIF periodic
      } // ENDELSE !InsideDomain
    } // ENDELSE

    if (DEBUG) 
      fprintf(stdout, "Walk: left grid: sent photon to grid %d (DeltaL = %d)\n", 
	      MoveToGrid->ID, DeltaLevel);
    break;

  case GravityUndefined:
  default:
    fprintf(stdout, "grid::WalkPhotonPackage: "
	    "GravityBoundaryType = RadiationBoundary undefined %"ISYM".\n",
	    GravityBoundaryType);
    ENZO_FAIL("");
  } // ENDSWITCH

  /* Error check */

#ifdef UNUSED
  if (MoveToGrid != NULL && InsideDomain)
    if (MoveToGrid->PointInGridNB(r) == FALSE)
      ENZO_FAIL("Photon not contained in MoveToGrid!");
#endif

  return SUCCESS;

}
