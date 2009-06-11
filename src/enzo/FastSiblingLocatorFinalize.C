/***********************************************************************
/
/  INITIALIZE THE CHAINING MESH FOR THE FAST SIBLING LOCATOR
/
/  written by: Greg Bryan
/  date:       October 2004
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
 
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh)
{
 
  int i, dim, size = 1;
  ChainingMeshLink *link, *oldlink;
 
  /* Calcualte number of mesh points in chaining mesh. */
 
  for (dim = 0; dim < Mesh->Rank; dim++)
    size *= Mesh->Dimension[dim];
 
  /* Delete linked list chains. */
 
  for (i = 0; i < size; i++) {
    link = Mesh->HeadOfChain[i];
    while (link != NULL) {
      oldlink = link;
      link = oldlink->NextLink;
      delete oldlink;
    }
  }
 
  /* Delete linked-list for grids outside the mesh. */
 
  link = Mesh->OutsideMesh;
  while (link != NULL) {
    oldlink = link;
      link = oldlink->NextLink;
    delete oldlink;
  }
 
  /* Delete the mesh itself. */
 
  delete [] Mesh->HeadOfChain;
 
  return SUCCESS;
}
