/***********************************************************************
/
/  GRID CLASS (CREATE AND FILL IN BOUNDARY FLUX STRUCTURE FOR THIS GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// Clear the boundary fluxes (allocate first if necesary).
//  Note that we assume the left and right flux boundaries are the same size.
//

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
 
void grid::PrepareBoundaryFluxes()
{
 
  int dim, i, j, field;
 
  /* If unallocated, then create BoundaryFluxes. */
 
  if (BoundaryFluxes == NULL)
    BoundaryFluxes = new fluxes;
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Set up the BoundaryFluxes descriptors (for each dimension [dim], there
       are two faces (left, right) that are each described by two vectors
       with indicies [j].  We set all the starting vectors to the left-most
       corner and the ending vectors to the right-most corner and then
       correct. */
 
    for (j = 0; j < GridRank; j++) {
      BoundaryFluxes->LeftFluxStartGlobalIndex[j][dim] =
	nlongint(( GridLeftEdge[dim] - DomainLeftEdge[dim]) /
		 CellWidth[dim][0]);
      BoundaryFluxes->LeftFluxEndGlobalIndex[j][dim] = max(
	 nlongint((GridRightEdge[dim] - DomainLeftEdge[dim]) /
		  CellWidth[dim][0])-1,
	 BoundaryFluxes->LeftFluxStartGlobalIndex[j][dim]);
      BoundaryFluxes->RightFluxStartGlobalIndex[j][dim] =
	  nlongint(( GridLeftEdge[dim] - DomainLeftEdge[dim]) /
		   CellWidth[dim][0]);
      BoundaryFluxes->RightFluxEndGlobalIndex[j][dim] = max(
	 nlongint((GridRightEdge[dim] - DomainLeftEdge[dim]) /
		  CellWidth[dim][0])-1,
	 BoundaryFluxes->RightFluxStartGlobalIndex[j][dim]);
    }
 
    /* correct so vectors point to correct position */
 
    BoundaryFluxes->LeftFluxEndGlobalIndex[dim][dim] =
      BoundaryFluxes->LeftFluxStartGlobalIndex[dim][dim];
    BoundaryFluxes->RightFluxStartGlobalIndex[dim][dim] =
      BoundaryFluxes->RightFluxEndGlobalIndex[dim][dim];
 
  } // end: loop over dims
 
  /* set Flux pointers to NULL */
 
  for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++)
    for (i = 0; i < MAX_DIMENSION; i++) {
      BoundaryFluxes->LeftFluxes[field][i]  = NULL;
      BoundaryFluxes->RightFluxes[field][i] = NULL;
    }
 
}
