/***********************************************************************
/
/  INITIALIZE ROOT GRID CHAINING MESH FOR THE FAST SIBLING LOCATOR
/
/  written by: Robert Harkness
/  date:       February, 2008
/  modified1:
/
/  PURPOSE:
/    In order to speed up the sibling search, we have implemented a chaining
/    mesh approach in order to localize the mesh and minimize the number
/    of grid-grid overlap checks.  This routine initializes the chaining mesh.
/    This parameters here could surely be tuned.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[])
{
 
  /* Fill out structure. */
 
  int dim, i, size = 1;
  for (dim = 0; dim < Rank; dim++) {
    Mesh->Rank = Rank;
    //    Mesh->Dimension[dim] = 64;  // this should be tuned
    //    Mesh->Dimension[dim] = TopGridDims[dim] / 4;
    Mesh->Dimension[dim] = 16;
 
    /* Set the chaining mesh size to be the same as the entire domain. */
 
    Mesh->LeftEdge[dim] = DomainLeftEdge[dim];
    Mesh->RightEdge[dim] = DomainRightEdge[dim];
    Mesh->BoundaryType[dim] = periodic;
 
    /* The other option is to just use the AMR refine region. */
 
    //    Mesh->LeftEdge[dim] = RefineRegionLeftEdge[dim];
    //    Mesh->RightEdge[dim] = RefineRegionRightEdge[dim];
    //    Mesh->BoundaryType[dim] = outflow; // i.e. not periodic
 
    Mesh->CellSize[dim] = (Mesh->RightEdge[dim] - Mesh->LeftEdge[dim])/
                          Mesh->Dimension[dim];
    size *= Mesh->Dimension[dim];
  }
 
  /* Generate empty head-of-chain stubs for each point in mesh. */
 
  Mesh->HeadOfChain = new ChainingMeshLink *[size];
  for (i = 0; i < size; i++)
    Mesh->HeadOfChain[i] = NULL;
  Mesh->OutsideMesh = NULL;
 
  return SUCCESS;
}
