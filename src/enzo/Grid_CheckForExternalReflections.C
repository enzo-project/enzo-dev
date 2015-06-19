/***********************************************************************
/
/  GRID CLASS (COPY ZONES ALONG EXTERNAL BOUNDARY, IF REFLECTING BCs)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Alexei Kritsuk, September 2005
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

// This routine checks whether this grid touches external boundary(ies)
// and copies zones (if boundary conditions are reflecting).

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

int grid::CheckForExternalReflections(boundary_type LeftFaceBoundaryCondition[],
				      boundary_type RightFaceBoundaryCondition[])
{

  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* declarations */

  int   i, j, k, dim, field, GhostZones = NumberOfGhostZones;
  float *index, Sign;
  int GridLeftOffset[MAX_DIMENSION], GridRightOffset[MAX_DIMENSION];


  /*
  if (debug)
    printf("Check for external reflections.\n");
  */

  /* This unfortunate bit of logic is to make sure we should be
     applying reflecting bc's in this direction. */

  for (dim = 0; dim < GridRank; dim++) {

    /* Left faces */

    GridLeftOffset[dim]  = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			         CellWidth[dim][0]);

    if (LeftFaceBoundaryCondition[dim] == reflecting &&
	GridLeftOffset[dim] < GhostZones)

      for (field = 0; field < NumberOfBaryonFields; field++) {
	Sign = 1;

	switch(dim) {

	/* copy x-left */

	case 0:
	  if (FieldType[field] == Velocity1) Sign = -1;
	  for (i = 0; i < GridStartIndex[0] - GridLeftOffset[0]; i++)
	    for (j = 0; j < GridDimension[1]; j++)
	      for (k = 0; k < GridDimension[2]; k++) {
		index = BaryonField[field] + i + 
		        j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
		*index = Sign*(*(index + (2*GridStartIndex[0] - 1 - 2*i)));
	      }
	  break;

	/* copy y-left */

	case 1:
	  if (FieldType[field] == Velocity2) Sign = -1;
	  for (j = 0; j < GridStartIndex[1] - GridLeftOffset[1]; j++)
	    for (i = 0; i < GridDimension[0]; i++)
	      for (k = 0; k < GridDimension[2]; k++) {
		index  = BaryonField[field] + i + 
		         j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
		*index = Sign*(*(index + (2*GridStartIndex[1] - 1 - 2*j)*GridDimension[0]));
	      }
	  break;

	/* copy z-left */

	case 2:      
	  if (FieldType[field] == Velocity3) Sign = -1;
	  for (k = 0; k < GridStartIndex[2] - GridLeftOffset[2]; k++)
	    for (i = 0; i < GridDimension[0]; i++)
	      for (j = 0; j < GridDimension[1]; j++) {
		index  = BaryonField[field] + i + 
		         j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
		*index = Sign*(*(index + (2*GridStartIndex[2]-1 - 2*k)*GridDimension[0]*
				 GridDimension[1]));
	  }	  
	  break;

	default:
	  fprintf(stderr, "G_CFER: Strange dimensionality...\n");
	  break;

	} /* switch */

      } /* loop over fields */

    /* Right faces */

    GridRightOffset[dim] = nint((DomainRightEdge[dim] - GridRightEdge[dim])/
			         CellWidth[dim][0]);

    if (RightFaceBoundaryCondition[dim] == reflecting &&
	GridRightOffset[dim] < GhostZones)

      for (field = 0; field < NumberOfBaryonFields; field++) {
	Sign = 1;

	switch(dim) {

	/* copy x-right */

	case 0:
	  if (FieldType[field] == Velocity1) Sign = -1;
	  for (i = 0; i < GridDimension[0]-GridEndIndex[0]-1-GridRightOffset[0]; i++)
	    for (j = 0; j < GridDimension[1]; j++)
	      for (k = 0; k < GridDimension[2]; k++) {
		index = BaryonField[field] + GridDimension[0] - 1 - i
		        + j*GridDimension[0]
		        + k*GridDimension[1]*GridDimension[0];
		*index = Sign*(*(index - 2 * (GridDimension[0] - GridEndIndex[0] - i) + 1));
	      }
	  break;

	/* copy y-right */

	case 1:
	  if (FieldType[field] == Velocity2) Sign = -1;
	  for (j = 0; j < GridDimension[1]-GridEndIndex[1]-1-GridRightOffset[1]; j++)
	    for (i = 0; i < GridDimension[0]; i++)
	      for (k = 0; k < GridDimension[2]; k++) {
		index = BaryonField[field] + i
		  + (GridDimension[1] - 1 - j)*GridDimension[0]
		  + k*GridDimension[1]*GridDimension[0];
		*index = Sign*(
		  *(index - (
		    2 * (GridDimension[1] - GridEndIndex[1] - j) + 1
		  ) * GridDimension[0])
		);
	      }
	  break;

	/* copy z-right */

	case 2:      
	  if (FieldType[field] == Velocity3) Sign = -1;
	  for (k = 0; k < GridDimension[2]-GridEndIndex[2]-1-GridRightOffset[2]; k++)
	    for (i = 0; i < GridDimension[0]; i++)
	      for (j = 0; j < GridDimension[1]; j++) {
		index = BaryonField[field] + i
		  + j*GridDimension[0]
		  + (
		    GridDimension[2] - 1 - k
		  ) * GridDimension[1] * GridDimension[0];
		*index = Sign*(
		  *(index - (
		    2 * (GridDimension[2] - GridEndIndex[2] - k) + 1
		  ) * GridDimension[0] * GridDimension[1])
		);
	      }	  
	  break;

	default:
	  fprintf(stderr, "G_CFER: Strange dimensionality...\n");
	  break;

	} /* switch */

    } /* loop over fields */

  } /* loop over dims */

  return SUCCESS;

}
