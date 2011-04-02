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
 
 
int grid::CheckForPossibleOverlap(grid *OtherGrid,
			  boundary_type LeftFaceBoundaryCondition[],
			  boundary_type RightFaceBoundaryCondition[])
{
 
  /* declarations */
 
  int i, j, k;
  FLOAT EdgeOffset[MAX_DIMENSION] = {0,0,0};
 
  /* Check all 26 neighbours. */
 
  int FullPeriod = TRUE;
 
  /* Always overlap with self. */

  if (this == OtherGrid)
    return TRUE;
 
  /* Check for overlap & return TRUE if there is some. */
 
  if (CheckForPossibleOverlapHelper(OtherGrid, EdgeOffset) == TRUE)
    return TRUE;
 
  /* For periodic boundary conditions, do some extra checks.  This insures
     that grids which overlap along periodic boundaries are handled correctly;
     we do this by actually moving the grid to it's periodic location in the
     next (periodic) domain over.
     (Here we use EdgeOffset to tell Grid_CopyZonesFromGrid that we have
     moved the base location of the grid). */
 

  FLOAT Lx, Ly, ShearingOffset;

  if (ShearingBoundaryDirection>-1) { // For shearing box we have another offset in the y direction
    Lx = (DomainRightEdge[ShearingBoundaryDirection]-DomainLeftEdge[ShearingBoundaryDirection]);
    Ly = (DomainRightEdge[ShearingVelocityDirection]-DomainLeftEdge[ShearingVelocityDirection]);
    ShearingOffset = AngularVelocity*VelocityGradient*Time*Lx;
    while (ShearingOffset >= Ly) {
      ShearingOffset -= Ly;
    }  
  }


  int kdim = (GridRank > 2) ? 1 : 0;
  int jdim = (GridRank > 1) ? 1 : 0;
  for (k = -kdim; k <= +kdim; k++) {
    EdgeOffset[2] = FLOAT(k)*(DomainRightEdge[2] - DomainLeftEdge[2]);
    for (j = -jdim; j <= +jdim; j++) {
      EdgeOffset[1] = FLOAT(j)*(DomainRightEdge[1] - DomainLeftEdge[1]);
      for (i = -1; i <= +1; i++) {
	EdgeOffset[0] = FLOAT(i)*(DomainRightEdge[0] - DomainLeftEdge[0]);
 
	/* This unfortunate bit of logic is to make sure we should be
	   applying periodic bc's in this direction. */


	if ((i != +1 || ((LeftFaceBoundaryCondition[0] == periodic || LeftFaceBoundaryCondition[0] == shearing) &&
			 (CellLeftEdge[0][0] < DomainLeftEdge[0] || ShearingVelocityDirection==0 ))    ) &&
	    (i != -1 || ((RightFaceBoundaryCondition[0] == periodic || RightFaceBoundaryCondition[0] == shearing) &&
			 (CellLeftEdge[0][GridDimension[0]-1] >
			 DomainRightEdge[0] ||  ShearingVelocityDirection==0 ))                        ) &&
	    (j != +1 || ((LeftFaceBoundaryCondition[1] == periodic || LeftFaceBoundaryCondition[1] == shearing) &&
			 (CellLeftEdge[1][0] < DomainLeftEdge[1] || ShearingVelocityDirection==1 ))    ) &&
	    (j != -1 || ((RightFaceBoundaryCondition[1] == periodic || RightFaceBoundaryCondition[1] == shearing) &&
			 (CellLeftEdge[1][GridDimension[1]-1] >
			 DomainRightEdge[1]  || ShearingVelocityDirection==1 ))                        ) &&
	    (k != +1 || ((LeftFaceBoundaryCondition[2] == periodic || LeftFaceBoundaryCondition[2] == shearing) &&
			 (CellLeftEdge[2][0] < DomainLeftEdge[2]  || ShearingVelocityDirection==2))    ) &&
	    (k != -1 || ((RightFaceBoundaryCondition[2] == periodic || RightFaceBoundaryCondition[2] == shearing) &&
			 (CellLeftEdge[2][GridDimension[2]-1] >
			 DomainRightEdge[2])  || ShearingVelocityDirection==2 )  )   ){
 
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
	


	  /* Full periodic case (26 checks). */
 
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
	    if (this->CheckForPossibleOverlapHelper(OtherGrid, EdgeOffset)
		== TRUE)
	      return TRUE;
	  }

	  EdgeOffset[2] = FLOAT(k)*(DomainRightEdge[2] - DomainLeftEdge[2]);
	  EdgeOffset[1] = FLOAT(j)*(DomainRightEdge[1] - DomainLeftEdge[1]);
	  EdgeOffset[0] = FLOAT(i)*(DomainRightEdge[0] - DomainLeftEdge[0]);
	 
	} // end: if (periodic bc's)


 
 
      } // end: loop of i
    } // end: loop of j
  } // end: loop of k
 
  /* If we get to this pointer, there is no overlap. */
 
  return FALSE;
 
}
 
int grid::CheckForPossibleOverlapHelper(grid *OtherGrid,
					FLOAT EdgeOffset[MAX_DIMENSION])
{
 
  /* declarations */
 
  int dim;
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
 
  /* Initialize gravitating mass field parameters. */

  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED)
    this->InitializeGravitatingMassField(RefineBy);
  if (OtherGrid->GravitatingMassFieldCellSize == FLOAT_UNDEFINED)
    OtherGrid->InitializeGravitatingMassField(RefineBy);
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    Left[dim] = max(
          GravitatingMassFieldLeftEdge[dim] + EdgeOffset[dim],
	  OtherGrid->GravitatingMassFieldLeftEdge[dim]);
 
    Right[dim] = min(
	  GravitatingMassFieldLeftEdge[dim] +
	  GravitatingMassFieldCellSize * GravitatingMassFieldDimension[dim] +
	  EdgeOffset[dim],
	  OtherGrid->GravitatingMassFieldLeftEdge[dim] +
	  OtherGrid->GravitatingMassFieldCellSize *
	  OtherGrid->GravitatingMassFieldDimension[dim]);
 
    if (Left[dim]+0.5*GravitatingMassFieldCellSize >= Right[dim])
      return FAIL;
  }
 
  return TRUE;
}
