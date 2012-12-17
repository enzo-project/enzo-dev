/***********************************************************************
/
/  GRID CLASS (Checks other grids for shared external faces)
/              Like this :
 
|--------|
|        |---|
|        |   |
|        |---|
|        |
|--------|
 
/
/  written by: Greg Bryan
/  date:       May, 1995
/  Re-Purposed for faces: June 2005, D Collins
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine checks for shared surfaces.
 
#include <math.h>
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
 
 
int grid::CheckForSharedFace(grid *OtherGrid,
			  boundary_type LeftFaceBoundaryCondition[],
			  boundary_type RightFaceBoundaryCondition[])
{
 
  /* declarations */
 
  int i, j, k, dim;
  FLOAT EdgeOffset[MAX_DIMENSION] = {0,0,0};
 
  /* Check all 26 neighbours. */
 
  int FullPeriod = TRUE;
 
  /* Always overlap with self. */
 
  if (this == OtherGrid)
    return TRUE;
 
  /* Check for overlap & return TRUE if there is some. */
 
  if (CheckForSharedFaceHelper(OtherGrid, EdgeOffset) == TRUE)
    return TRUE;
 
  /* For periodic boundary conditions, do some extra checks.  This insures
     that grids which overlap along periodic boundaries are handled correctly;
     we do this by actually moving the grid to it's periodic location in the
     next (periodic) domain over.
     (Here we use EdgeOffset to tell Grid_CopyZonesFromGrid that we have
     moved the base location of the grid). */
 
  int kdim = (GridRank > 2) ? 1 : 0;
  int jdim = (GridRank > 1) ? 1 : 0;
 

  FLOAT Lx, Ly, ShearingOffset;

  if (ShearingBoundaryDirection>-1) { // For shearing box we have another offset in the y direction
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
    BoundaryCheck[2*dim] = false;
    BoundaryCheck[2*dim+1] = false;
  }

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

  for (k = -kdim; k <= +kdim; k++) {
    EdgeOffset[2] = FLOAT(k) * DomainWidth[2];
    zcheck = (k != +1 || BoundaryCheck[4]) && (k != -1 || BoundaryCheck[5]);
    for (j = -jdim; j <= +jdim; j++) {
      EdgeOffset[1] = FLOAT(j) * DomainWidth[1];
      ycheck = (j != +1 || BoundaryCheck[2]) && (j != -1 || BoundaryCheck[3]);
      for (i = -1; i <= +1; i++) {
	EdgeOffset[0] = FLOAT(i) * DomainWidth[0];
 
	/* This unfortunate bit of logic (JHW, Dec09: moved outside
	   the loop) is to make sure we should be applying periodic
	   bc's in this direction. */
 
	if ((i != +1 || BoundaryCheck[0]) &&
	    (i != -1 || BoundaryCheck[1]) &&
	    ycheck && zcheck) {

// 	if ((i != +1 || ((LeftFaceBoundaryCondition[0] == periodic || LeftFaceBoundaryCondition[0] == shearing) &&
// 			 (CellLeftEdge[0][0] < DomainLeftEdge[0] ))    ) &&
// 	    (i != -1 || ((RightFaceBoundaryCondition[0] == periodic || RightFaceBoundaryCondition[0] == shearing) &&
// 			 (CellLeftEdge[0][GridDimension[0]-1] >
// 			 DomainRightEdge[0] ))                        ) &&
// 	    (j != +1 || ((LeftFaceBoundaryCondition[1] == periodic || LeftFaceBoundaryCondition[1] == shearing) &&
// 			 (CellLeftEdge[1][0] < DomainLeftEdge[1]  ))    ) &&
// 	    (j != -1 || ((RightFaceBoundaryCondition[1] == periodic || RightFaceBoundaryCondition[1] == shearing) &&
// 			 (CellLeftEdge[1][GridDimension[1]-1] >
// 			 DomainRightEdge[1] ))                        ) &&
// 	    (k != +1 || ((LeftFaceBoundaryCondition[2] == periodic || LeftFaceBoundaryCondition[2] == shearing) &&
// 			 (CellLeftEdge[2][0] < DomainLeftEdge[2] ))    ) &&
// 	    (k != -1 || ((RightFaceBoundaryCondition[2] == periodic || RightFaceBoundaryCondition[2] == shearing) &&
// 			 (CellLeftEdge[2][GridDimension[2]-1] >
// 			 DomainRightEdge[2])  )  )   ){

 
	  /* Full periodic case (26 checks).
	     This ONLY checks the Periodic shifts.  (that's the i!=0 || ... crap) */
	


	  if ((GridRank > 2 || k == 0) &&
	      (GridRank > 1 || j == 0) &&
	      (i != 0 || j != 0 || k != 0)) {

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
 
	    if (this->CheckForSharedFaceHelper(OtherGrid, EdgeOffset)
		== TRUE)
	      return TRUE;
	  }
 
	  if (ShearingBoundaryDirection!=-1){
	    EdgeOffset[2] = FLOAT(k)*DomainWidth[2];
	    EdgeOffset[1] = FLOAT(j)*DomainWidth[1];
	    EdgeOffset[0] = FLOAT(i)*DomainWidth[0];
	  }

	} // end: if (periodic bc's)

 
 
      } // end: loop of i
    } // end: loop of j
  } // end: loop of k
 
  /* If we get to this pointer, there is no overlap. */
 
  return FALSE;
 
}
 
int grid::CheckForSharedFaceHelper(grid *OtherGrid,
					FLOAT EdgeOffset[MAX_DIMENSION])
{
  /* declarations */
 
 
  int dim, dim2, GotSide=0;
  FLOAT CellEpsilon[MAX_DIMENSION];
 
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++) {
    CellEpsilon[dim] = min( 0.5*CellWidth[dim][0],
			    0.5*OtherGrid->CellWidth[dim][0] );
    Left[dim] = max(
		    GridLeftEdge[dim] + EdgeOffset[dim],
		    OtherGrid->GridLeftEdge[dim]);
 
    Right[dim] = min(
		     GridRightEdge[dim] + EdgeOffset[dim],
		     OtherGrid->GridRightEdge[dim]);
 
    //If there's any space, fail
 
    if (Left[dim]-CellEpsilon[dim] >= Right[dim])
      return FALSE;
  }//dim
 
  //here goes the rest of the check
 
  for(dim=0;dim<GridRank;dim++){
    if( fabs(Left[dim]-Right[dim]) < CellEpsilon[dim] ){
 
      for(dim2=0;dim2<GridRank;dim2++)
	if(dim2 != dim){
	  if( fabs(Left[dim2]-Right[dim2]) < CellEpsilon[dim2] )
	    return FALSE;
	}//third dim loop
 
      return TRUE;
 
    }//if face matches
  }//dim again
 
    return FALSE;
}
 
