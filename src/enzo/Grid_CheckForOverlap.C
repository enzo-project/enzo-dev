/***********************************************************************
/
/  GRID CLASS (CHECK OTHER GRID FOR ANY OVERLAPPING ZONES)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine checks the other grid for any overlapping zones (including
//   periodic boundary conditions).  If any are found, CopyZonesFromGrid
//   is called.
 
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int grid::CheckForOverlap(grid *OtherGrid,
			  boundary_type LeftFaceBoundaryCondition[],
			  boundary_type RightFaceBoundaryCondition[],
			  int (grid::*CopyFunction)(grid *OtherGrid,
						    FLOAT EdgeOffset[]))
{
 
  // Return if this doesn't involve us
 
  if (this->CommunicationMethodShouldExit(OtherGrid))
    return SUCCESS;
 
  int i, j, k, dim;
  FLOAT EdgeOffset[MAX_DIMENSION] = {0.0,0.0,0.0};
 
  /* If the copy function is AddOverlappingParticleMassField, then
     apply to self, otherwise don't. */
  
#ifndef TRANSFER
  int DoSelf = (CopyFunction == &grid::AddOverlappingParticleMassField)?
    TRUE : FALSE;
#else 
  int DoSelf = (CopyFunction == &grid::AddOverlappingParticleMassField ||
		CopyFunction == &grid::SetSubgridMarkerFromSibling)?
    TRUE : FALSE;
#endif

  int FullPeriod = (CopyFunction == &grid::CopyPotentialField)?
    TRUE : FALSE;
 
  //  if (CopyFunction == &grid::CopyZonesFromGrid)
  //    FullPeriod = TRUE;
 
  //  Mod by RH - always check full periodic
 
   FullPeriod = TRUE;
 
  /*
    if (CopyFunction == &grid::AddOverlappingParticleMassField) {
    fprintf(stderr, "CFOverlap: AddOverlappingParticleMassField\n");
    }
 
    if (CopyFunction == &grid::CopyPotentialField) {
    fprintf(stderr, "CFOverlap: CopyPotentialField\n");
    }
 
    if (CopyFunction == &grid::CopyZonesFromGrid) {
    fprintf(stderr, "CFOverlap: CopyZonesFromGrid\n");
    }
  */
 
  // Check for overlap & do copy

   // TA HACK   This seems to get the correct GravitatingMassField .... ! 
   //           velocities look weird now ... more checking!
   //      DoSelf = FALSE;

  if (this != OtherGrid || DoSelf)
    if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
      ENZO_FAIL("Error in grid->*CopyFunction\n");
    }
 

  FLOAT Lx, Ly, ShearingOffset;

  if (ShearingBoundaryDirection!=-1) { // For shearing box we have another offset in the y direction
    Lx = (DomainRightEdge[ShearingBoundaryDirection]-DomainLeftEdge[ShearingBoundaryDirection]);
    Ly = (DomainRightEdge[ShearingVelocityDirection]-DomainLeftEdge[ShearingVelocityDirection]);
    ShearingOffset = AngularVelocity*VelocityGradient*Time*Lx;
    while (ShearingOffset >= Ly) {
      ShearingOffset -= Ly;
    }  
  }

  /* Pre-compute boundary checks for periodic bc's */

  bool BoundaryCheck[2*MAX_DIMENSION];
  bool ycheck, zcheck;
  FLOAT DomainWidth[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {

    BoundaryCheck[2*dim] = 
      ((LeftFaceBoundaryCondition[dim] == periodic || 
	LeftFaceBoundaryCondition[dim] == shearing) &&
       (CellLeftEdge[dim][0] < DomainLeftEdge[dim] || 
	ShearingVelocityDirection==dim ));

    BoundaryCheck[2*dim+1] = 
      ((RightFaceBoundaryCondition[dim] == periodic || 
	RightFaceBoundaryCondition[dim] == shearing) &&
       (CellLeftEdge[dim][GridDimension[dim]-1] > DomainRightEdge[dim] ||
	ShearingVelocityDirection==dim ));

    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    BoundaryCheck[2*dim] = TRUE;
    BoundaryCheck[2*dim+1] = TRUE;
    DomainWidth[dim] = 0.0;
  }

  /* For periodic boundary conditions, do some extra checks.  This insures
     that grids which overlap along periodic boundaries are handled correctly;
     we do this by actually moving the grid to it's periodic location in the
     next (periodic) domain over.
     (Here we use EdgeOffset to tell Grid_CopyZonesFromGrid that we have
     moved the base location of the grid). */
 

  //PrintToScreenBoundaries(BaryonField[2], "Vy");

  int kdim = (GridRank > 2) ? 1 : 0;
  int jdim = (GridRank > 1) ? 1 : 0;
  for (k = -kdim; k <= +kdim; k++) {
    EdgeOffset[2] = FLOAT(k) * DomainWidth[2];
    zcheck = (k != +1 || BoundaryCheck[4]) && (k != -1 || BoundaryCheck[5]);
    for (j = -jdim; j <= +jdim; j++) {
      EdgeOffset[1] = FLOAT(j) * DomainWidth[1];
      ycheck = (j != +1 || BoundaryCheck[2]) && (j != -1 || BoundaryCheck[3]);
      for (i = -1; i <= +1; i++) {
	EdgeOffset[0] = FLOAT(i) * DomainWidth[0];
 
	/* This unfortunate bit of logic (JHW Dec09: moved outside the
	   loop) is to make sure we should be applying periodic bc's
	   in this direction. */
 
	if ((i != +1 || BoundaryCheck[0]) &&
	    (i != -1 || BoundaryCheck[1]) &&
	    ycheck && zcheck) {

	  // Full periodic case (26 checks)
	    
	  if ((GridRank > 2 || k == 0) && 
	      (GridRank > 1 || j == 0) &&
	      (i != 0 || j != 0 || k != 0) && 
	      (FullPeriod == TRUE || ShearingBoundaryDirection!=-1)) {
	    	    
		if (ShearingBoundaryDirection!=-1){
		  if ((i== +1 && LeftFaceBoundaryCondition[0] == shearing) ||
		      (j== +1 && LeftFaceBoundaryCondition[1] == shearing) ||
		      (k== +1 && LeftFaceBoundaryCondition[2] == shearing)){
		    EdgeOffset[ShearingVelocityDirection] -= ShearingOffset;
		  }
		  if ((i== -1 && RightFaceBoundaryCondition[0] == shearing) ||
		      (j== -1 && RightFaceBoundaryCondition[1] == shearing) ||
		      (k== -1 && RightFaceBoundaryCondition[2] == shearing)){
		    EdgeOffset[ShearingVelocityDirection] += ShearingOffset;
		  }
		}

	    
		if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
		  printf("Error in grid->*CopyFunction (2)\n");
		  ENZO_FAIL("CopyFunctionFail(2)");
		}
	  }

	  // partial periodic case (6 checks)
	    
	  if ((GridRank > 2 || k == 0) && (GridRank > 1 || j == 0) &&
	      (FullPeriod==FALSE && ShearingBoundaryDirection==-1) && (ABS(i)+ABS(j)+ABS(k) == 1)) {

	    
	    if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
	      ENZO_FAIL("Error in grid->*CopyFunction (3)\n");
	    }
	  }

	  if (ShearingBoundaryDirection != -1) {
	    EdgeOffset[2] = FLOAT(k)*DomainWidth[2];
	    EdgeOffset[1] = FLOAT(j)*DomainWidth[1];
	    EdgeOffset[0] = FLOAT(i)*DomainWidth[0];
	  }

	} // end: if (periodic bc's)

 
      } // end: loop of i

    } // end: loop of j
  } // end: loop of k


  return SUCCESS;
 
}




 
