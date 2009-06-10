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
 
  int i, j, k;
  FLOAT EdgeOffset[MAX_DIMENSION] = {0,0,0};
 
  /* If the copy function is AddOverlappingParticleMassField, then
     apply to self, otherwise don't. */
 
  int DoSelf = (CopyFunction == &grid::AddOverlappingParticleMassField)?
		TRUE : FALSE;
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
 
  if (this != OtherGrid || DoSelf)
    if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
      fprintf(stderr, "Error in grid->*CopyFunction\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
 
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
 
	  // Full periodic case (26 checks)
 
	  if ((GridRank > 2 || k == 0) && (GridRank > 1 || j == 0) &&
	      (i != 0 || j != 0 || k != 0) && FullPeriod == TRUE) {
	    if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
	      fprintf(stderr, "Error in grid->*CopyFunction (2)\n");
	      ENZO_FAIL("Error in: "__FILE__);
	    }
	  }
 
	  // partial periodic case (6 checks)
	
	  if ((GridRank > 2 || k == 0) && (GridRank > 1 || j == 0) &&
	      FullPeriod == FALSE && (ABS(i)+ABS(j)+ABS(k) == 1)) {
	    if ((this->*CopyFunction)(OtherGrid, EdgeOffset) == FAIL) {
	      fprintf(stderr, "Error in grid->*CopyFunction (3)\n");
	      ENZO_FAIL("Error in: "__FILE__);
	    }
	  }
 
	} // end: if (periodic bc's)
 
      } // end: loop of i
    } // end: loop of j
  } // end: loop of k
 
  return SUCCESS;
 
}
