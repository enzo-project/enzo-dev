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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
#ifdef FLUX_FIX
 
int grid::CheckForSharedFace(grid *OtherGrid,
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
 
  for (k = -kdim; k <= +kdim; k++) {
    EdgeOffset[2] = FLOAT(k)*(DomainRightEdge[2] - DomainLeftEdge[2]);
    for (j = -jdim; j <= +jdim; j++) {
      EdgeOffset[1] = FLOAT(j)*(DomainRightEdge[1] - DomainLeftEdge[1]);
      for (i = -1; i <= +1; i++) {
	EdgeOffset[0] = FLOAT(i)*(DomainRightEdge[0] - DomainLeftEdge[0]);
 
	/* This unfortunate bit of logic is to make sure we should be
	   applying periodic bc's in this direction. */
 
	if ((i != +1 || (LeftFaceBoundaryCondition[0] == periodic &&
			 CellLeftEdge[0][0] < DomainLeftEdge[0])    ) &&
	    (i != -1 || (RightFaceBoundaryCondition[0] == periodic &&
			 CellLeftEdge[0][GridDimension[0]-1] >
			 DomainRightEdge[0])                        ) &&
	    (j != +1 || (LeftFaceBoundaryCondition[1] == periodic &&
			 CellLeftEdge[1][0] < DomainLeftEdge[1])    ) &&
	    (j != -1 || (RightFaceBoundaryCondition[1] == periodic &&
			 CellLeftEdge[1][GridDimension[1]-1] >
			 DomainRightEdge[1])                        ) &&
	    (k != +1 || (LeftFaceBoundaryCondition[2] == periodic &&
			 CellLeftEdge[2][0] < DomainLeftEdge[2])    ) &&
	    (k != -1 || (RightFaceBoundaryCondition[2] == periodic &&
			 CellLeftEdge[2][GridDimension[2]-1] >
			 DomainRightEdge[2])                        )   ) {
 
	  /* Full periodic case (26 checks).
	     This ONLY checks the Periodic shifts.  (that's the i!=0 || ... crap) */
	
	  if ((GridRank > 2 || k == 0) &&
	      (GridRank > 1 || j == 0) &&
	      (i != 0 || j != 0 || k != 0)) {
 
	    if (this->CheckForSharedFaceHelper(OtherGrid, EdgeOffset)
		== TRUE)
	      return TRUE;
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
      return FAIL;
  }//dim
 
  //here goes the rest of the check
 
  for(dim=0;dim<GridRank;dim++){
    if( fabs(Left[dim]-Right[dim]) < CellEpsilon[dim] ){
 
      for(dim2=0;dim2<GridRank;dim2++)
	if(dim2 != dim){
	  if( fabs(Left[dim2]-Right[dim2]) < CellEpsilon[dim2] )
	    return FAIL;
	}//third dim loop
 
      return TRUE;
 
    }//if face matches
  }//dim again
 
    return FAIL;
}
 
#endif
