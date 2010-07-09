/***********************************************************************
/
/  GRID CLASS (FIND THE LIST OF POSSIBLE SIBLINGS FROM CHAINING MESH)
/
/  written by: Greg Bryan
/  date:       October, 2004
/  modified1:
/
/  PURPOSE: This routine looks up the list of possible siblings from
/     the chaining mesh, ignoring cases for which both grids are not
/     on this processor.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
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
 

 
int grid::FastSiblingLocatorFindSiblings(ChainingMeshStructure *Mesh,
			  SiblingGridList *list,
			  boundary_type LeftBoundaryCondition[],
			  boundary_type RightBoundaryCondition[])
{
 
  /* declarations */
 
  int index, n, i1,j1,k1, i, j, k, dim, istart[] = {0,0,0}, iend[] = {0,0,0};
  FLOAT Left, Right;
  ChainingMeshLink *current_link;
  static grid *TempList[MAX_NUMBER_OF_SUBGRIDS];

 
  /* Initialize gravitating mass field parameters. */
 
  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED)
    this->InitializeGravitatingMassField(RefineBy);
 
  /* Initialize the sibling grid list. */
 
  list->NumberOfSiblings = 0;
  list->GridList = NULL;
 
  /* Compute start and end indices of. */
 
  int OutsideMesh = FALSE, InsideMesh = TRUE, AlreadyPresent;
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Compute Left and Right edge of this grid.
       Use the GravitatingMassFieldLeftEdge and RightEdge because
       it is the largest region (and equal to CellLeft and Right if gravity
       is off). */
 
    Left = GravitatingMassFieldLeftEdge[dim];
    Right = GravitatingMassFieldLeftEdge[dim] +
	  GravitatingMassFieldCellSize * GravitatingMassFieldDimension[dim];
 
    /* If the boundary conditions are not periodic and this grid extends
       outside the chaining mesh then truncate it at the boundary and
       set the outside mesh flag. */
 
    if (Left < Mesh->LeftEdge[dim] && Mesh->BoundaryType[dim] != periodic) {
      Left = Mesh->LeftEdge[dim];
      OutsideMesh = TRUE;
    }
 
    if (Right > Mesh->RightEdge[dim] && Mesh->BoundaryType[dim] != periodic) {
      Right = Mesh->RightEdge[dim];
      OutsideMesh = TRUE;
    }
 
    /* If the grid is outside the chaining mesh entirely, then skip rest. */
 
    if (Right < Mesh->LeftEdge[dim] || Left > Mesh->RightEdge[dim]) {
      InsideMesh = FALSE;
      break;
    }
 
    /* Convert Left and Right into a (floating point) index and add the
       dimension size to it in case the Left or Right are negative
       (we take the mod of this later so it has no effect otherwise). */
 
    Left  = (Left  - Mesh->LeftEdge[dim])/Mesh->CellSize[dim] +
            Mesh->Dimension[dim];
    Right = (Right - Mesh->LeftEdge[dim])/Mesh->CellSize[dim] +
            Mesh->Dimension[dim];
 
    /* Convert into integer. */
 
    istart[dim] = nint(Left);
    iend[dim]   = nint(Right);
 
  }
 
  /* If this grid lies inside the chaining mesh, then loop over all
     mesh points that overlap with this grid and look up the grids. */
 
  if (InsideMesh) {
    for (k = istart[2]; k <= iend[2]; k++) {
      k1 = (GridRank > 2) ? k % Mesh->Dimension[2] : k;
      for (j = istart[1]; j <= iend[1]; j++) {
	j1 = (GridRank > 1) ? j % Mesh->Dimension[1] : j;
	index = (k1*Mesh->Dimension[1] + j1)*Mesh->Dimension[0];
	for (i = istart[0]; i <= iend[0]; i++) {
	  i1 = i % Mesh->Dimension[0];
	  current_link = Mesh->HeadOfChain[index+i1];
	  while (current_link != NULL) {
 
	    AlreadyPresent = FALSE;
            if ( list->NumberOfSiblings > MAX_NUMBER_OF_SUBGRIDS ) {
              ENZO_VFAIL("DC no of sibs > MAX_NUMBER_OF_SUBGRIDS %"ISYM"\n", list->NumberOfSiblings)
            }
	    for (n = 0; n < list->NumberOfSiblings; n++)
	      if (current_link->GridData == TempList[n])
		AlreadyPresent = TRUE;
 
	    /* If either grid is on this processor and there is some
	       possible overlap then add it to the temporary list. */
 
	    if (!AlreadyPresent &&
		(this->ProcessorNumber == MyProcessorNumber ||
		 current_link->GridData->ProcessorNumber == MyProcessorNumber)) {
	      if (this->CheckForPossibleOverlap(current_link->GridData,
					 LeftBoundaryCondition,
					 RightBoundaryCondition) == TRUE) {
                if ( list->NumberOfSiblings > MAX_NUMBER_OF_SUBGRIDS ) {
                  ENZO_VFAIL("DC2 no of sibs > MAX_NUMBER_OF_SUBGRIDS %"ISYM"\n", list->NumberOfSiblings)
                }
		TempList[list->NumberOfSiblings++] = current_link->GridData;
	      } 
	    }
 
	    /* Next Link in list. */
 
	    current_link = current_link->NextLink;
 
	  } // end: while (current_link != NULL)
 
	} // end: loop over i
      } // end: loop over j
    } // end: loop over k
 
  } // end: if (InsideMesh)
 
  /* If this grid lies outside the chaining mesh (either partially or
     entirely), then add this grid to that linked list. */
 
  if (OutsideMesh) {
    current_link = Mesh->OutsideMesh;
    while (current_link != NULL) {
 
      /* If either grid is on this processor and there is some
	 possible overlap then add it to the temporary list. */
 
      if (this->ProcessorNumber == MyProcessorNumber ||
	  current_link->GridData->ProcessorNumber == MyProcessorNumber) {
	if (this->CheckForPossibleOverlap(current_link->GridData,
					  LeftBoundaryCondition,
					  RightBoundaryCondition) == TRUE) {
          if ( list->NumberOfSiblings > MAX_NUMBER_OF_SUBGRIDS ) {
            ENZO_VFAIL("DC3 no of sibs > MAX_NUMBER_OF_SUBGRIDS %"ISYM"\n", list->NumberOfSiblings)
          }

	  TempList[list->NumberOfSiblings++] = current_link->GridData;
	}
      }
 
      /* Next Link in list. */
 
      current_link = current_link->NextLink;
 
    } // end: while (current_link != NULL)
 
  }
 
  /* Create an array and copy over list of grids. */
 
  if (list->NumberOfSiblings > 0) {
    list->GridList = new grid *[list->NumberOfSiblings];


    for (i = 0; i < list->NumberOfSiblings; i++) {
      if ( i > MAX_NUMBER_OF_SUBGRIDS ) {

        fprintf(stderr, "DC4 no of sibs > MAX_NUMBER_OF_SUBGRIDS %"ISYM"\n", i);
      }
      list->GridList[i] = TempList[i];
    }
  }
 
  return SUCCESS;
}
