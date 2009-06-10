/***********************************************************************
/
/  GRID CLASS (ADD A GRID INTO THE CHAINING MESH)
/
/  written by: Greg Bryan
/  date:       October, 2004
/  modified1:
/
/  PURPOSE: This routine adds links into the part of the chaining mesh
/     which the grid overlays.
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

 
int grid::FastSiblingLocatorAddGrid(ChainingMeshStructure *Mesh)
{
 
  /* declarations */
 
  int index, i, j, k, i1, j1, k1, dim, istart[] = {0,0,0}, iend[] = {0,0,0};
  FLOAT Left, Right;
  ChainingMeshLink *new_link;
 
  /* Initialize gravitating mass field parameters. */
 
  if (GravitatingMassFieldCellSize == FLOAT_UNDEFINED)
    this->InitializeGravitatingMassField(RefineBy);
 
  /* Compute start and end indices of. */
 
  int OutsideMesh = FALSE, InsideMesh = TRUE;
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
     mesh points that overlap with this grid and add a linked list pointing
     to this grid. */
 
  if (InsideMesh) {
    for (k = istart[2]; k <= iend[2]; k++) {
      k1 = (GridRank > 2) ? k % Mesh->Dimension[2] : k;
      for (j = istart[1]; j <= iend[1]; j++) {
	j1 = (GridRank > 1) ? j % Mesh->Dimension[1] : j;
	index = (k1*Mesh->Dimension[1] + j1)*Mesh->Dimension[0];
	for (i = istart[0]; i <= iend[0]; i++) {
	  i1 = i % Mesh->Dimension[0];
	  new_link = new ChainingMeshLink;
	  new_link->NextLink = Mesh->HeadOfChain[index+i1];
	  new_link->GridData = this;
	  Mesh->HeadOfChain[index+i1] = new_link;
	}
      }
    }
  }
 
  /* If this grid lies outside the chaining mesh (either partially or
     entirely), then add this grid to that linked list. */
 
  if (OutsideMesh) {
    new_link = new ChainingMeshLink;
    new_link->NextLink = Mesh->OutsideMesh;
    new_link->GridData = this;
    Mesh->OutsideMesh = new_link;
  }
 
  return SUCCESS;
}
