/***********************************************************************
/
/  SetSubgridMarker FUNCTION
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:  John Wise, April 2010 - sets ghost zones with sibling 
/              grids.  If no sibling, set it to the parent grid.
/
/  PURPOSE: Set the SubgridMarker field for all grids on finer 
/           levels than this one
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

#define MIN_LEVEL 0

int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);

int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level)
{

  if (!RadiativeTransfer)
    return SUCCESS;

  int i, grid2;
  LevelHierarchyEntry *Temp;
  HierarchyEntry *Subgrid;
  int NumberOfGrids;

  ChainingMeshStructure ChainingMesh;
  SiblingGridList SiblingList;

  for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++)  {

    /* Initialize and fill out the fast sibling chaining mesh.  Only
       for AMR grids. */

    if (i >= MIN_LEVEL) {
      FastSiblingLocatorInitialize(&ChainingMesh, MetaData.TopGridRank,
				   MetaData.TopGridDims);
      for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);
    }

    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel) {

      /* 1. First the grid marks itself */

      Temp->GridData->SetSubgridMarkerFromSubgrid(Temp->GridData);
      
      /* 2. Mark the parent in the ghost zones */

      if (i > 0)
	Temp->GridData->SetSubgridMarkerFromParent
	  (Temp->GridHierarchyEntry->ParentGrid->GridData);
      //Temp->GridData->SetSubgridMarkerFromParent(NULL);

      /* 3. Mark subgrids next */

      Subgrid = Temp->GridHierarchyEntry->NextGridNextLevel;
      while (Subgrid != NULL) {
	Temp->GridData->SetSubgridMarkerFromSubgrid(Subgrid->GridData);
	Subgrid = Subgrid->NextGridThisLevel;
      } // ENDWHILE Subgrid

      /* 4. Mark ghost zones with any siblings for level>0 because we
	 have to treat the domain boundaries differntly */
      
      /* Get a list of possible siblings from the chaining mesh */

      if (i >= MIN_LEVEL) {
	Temp->GridData->FastSiblingLocatorFindSiblings
	  (&ChainingMesh, &SiblingList, MetaData.LeftFaceBoundaryCondition,
	   MetaData.RightFaceBoundaryCondition);

	for (grid2 = 0; grid2 < SiblingList.NumberOfSiblings; grid2++)
	  if (Temp->GridData != SiblingList.GridList[grid2])
	    Temp->GridData->
	      CheckForOverlap(SiblingList.GridList[grid2],
			      MetaData.LeftFaceBoundaryCondition,
			      MetaData.RightFaceBoundaryCondition,
			      &grid::SetSubgridMarkerFromSibling);

	delete [] SiblingList.GridList;
      } // ENDIF i > 0

    } // ENDFOR grids

    if (i >= MIN_LEVEL)
      FastSiblingLocatorFinalize(&ChainingMesh);      

  } // ENDFOR levels

  return SUCCESS;
}
