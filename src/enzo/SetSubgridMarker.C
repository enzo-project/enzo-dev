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
#ifdef USE_MPI
#include "mpi.h"
#endif
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
#include "communication.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);

int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level,
		     int UpdateReplicatedGridsOnly)
{

  if (!RadiativeTransfer)
    return SUCCESS;

  int i, grid1, grid2, StartLevel;
  LevelHierarchyEntry *Temp;
  HierarchyEntry *Subgrid;

  ChainingMeshStructure ChainingMesh;
  SiblingGridList SiblingList;
  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int NumberOfGrids[MAX_DEPTH_OF_HIERARCHY];
  char *UpdateMask[MAX_DEPTH_OF_HIERARCHY];

  /* Generate grid arrays for lookups */

  for (i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++) {
    UpdateMask[i] = NULL;
    if (LevelArray[i] != NULL)
      NumberOfGrids[i] = GenerateGridArray(LevelArray, i, &Grids[i]);
    else
      NumberOfGrids[i] = 0;
  }

  // Initialize SubgridMarker only on root grid
  if (level == 0)
    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->SetSubgridMarkerFromSubgrid(NULL);

  StartLevel = max(level, 0);

  /* Generate update masks.  Only update the replicated grids (in the
     RT load balancing) if the argument UpdateReplicatedGridsOnly is
     true. */

  int temp_proc, ori_proc;

  if (UpdateReplicatedGridsOnly) {
    for (i = StartLevel; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {
      if (NumberOfGrids[i] > 0) {
	UpdateMask[i] = new char[NumberOfGrids[i]];
	for (grid1 = 0; grid1 < NumberOfGrids[i]; grid1++) {
	  temp_proc = Grids[i][grid1]->GridData->ReturnProcessorNumber();
	  ori_proc = Grids[i][grid1]->GridData->ReturnOriginalProcessorNumber();
	  UpdateMask[i][grid1] = (temp_proc != ori_proc) ? TRUE : FALSE;
	} // ENDFOR grids
      } // ENDIF grids>0
    } // ENDFOR level
  } else {
    for (i = StartLevel; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {
      if (NumberOfGrids[i] > 0) {
	UpdateMask[i] = new char[NumberOfGrids[i]];
	for (grid1 = 0; grid1 < NumberOfGrids[i]; grid1++)
	  UpdateMask[i][grid1] = TRUE;
      } // ENDIF grids>0
    } // ENDFOR level
  } // ENDELSE

  for (i = StartLevel; i < MAX_DEPTH_OF_HIERARCHY-1; i++)  {

    /* Initialize and fill out the fast sibling chaining mesh. */

    FastSiblingLocatorInitialize(&ChainingMesh, MetaData.TopGridRank,
				 MetaData.TopGridDims);
    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

    /* 1. First the grid marks itself */

    for (Temp = LevelArray[i], grid1 = 0; Temp; 
	 Temp = Temp->NextGridThisLevel, grid1++)
      if (UpdateMask[i][grid1]) {
	Temp->GridData->
	  CheckForOverlap(Temp->GridData,
			  MetaData.LeftFaceBoundaryCondition,
			  MetaData.RightFaceBoundaryCondition,
			  &grid::SetSubgridMarkerFromSibling);
      }
      //Temp->GridData->SetSubgridMarkerFromSubgrid(Temp->GridData);
      
    /* 2. Mark the parent in the ghost zones.  MOVED: after the grid
       loop to include SubgridMarker in parent's ghost zones. */

    /******************************************************************
      Some of the ghost zones might be covered by grids other than
      siblings and parents in the situation when the parent and child
      boundaries match.  Thus the ghost zones will be some other
      coarse grid.  Instead of searching through the grids again, this
      information is in the ghost zones of the parent's SubgridMarker,
      which may be on a different processor.  Retrieve these values.
    *******************************************************************/

    if (i > 0) {
      
      /* First stage: Post receives */

      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;
      for (Temp = LevelArray[i], grid1=0; Temp; 
	   Temp = Temp->NextGridThisLevel, grid1++) {
	if (UpdateMask[i][grid1])
	  Temp->GridData->SetSubgridMarkerFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData, i);
      }

      /* Second stage: Send data */

      CommunicationDirection = COMMUNICATION_SEND;
      for (Temp = LevelArray[i], grid1=0; Temp; 
	   Temp = Temp->NextGridThisLevel, grid1++) {
	if (UpdateMask[i][grid1])
	  Temp->GridData->SetSubgridMarkerFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData, i);
      }

      /* Third stage: Receive and process data */

      CommunicationReceiveHandler();

      /* In parallel, we need to convert the packed integers (stored
	 in the ghost zones of SubgridMarker) received during
	 communication into grid pointers if the parent is on a
	 different processor */
      
      for (Temp = LevelArray[i], grid1=0; Temp; 
	   Temp = Temp->NextGridThisLevel, grid1++)
	if (UpdateMask[i][grid1])
	  Temp->GridData->SubgridMarkerPostParallelGZ
	    (Temp->GridHierarchyEntry->ParentGrid->GridData, Grids,
	     NumberOfGrids);

    } // ENDIF level > 0

    /* 3. Mark subgrids next */

    for (Temp = LevelArray[i], grid1=0; Temp; 
	 Temp = Temp->NextGridThisLevel, grid1++)
      if (UpdateMask[i][grid1]) {
	Subgrid = Temp->GridHierarchyEntry->NextGridNextLevel;
	while (Subgrid != NULL) {
	  Temp->GridData->SetSubgridMarkerFromSubgrid(Subgrid->GridData);
	  Subgrid = Subgrid->NextGridThisLevel;
	} // ENDWHILE Subgrid
      } // ENDIF UpdateMask

    /* 4. Mark ghost zones with any siblings for level>0 because we
       have to treat the domain boundaries differntly */

    for (Temp = LevelArray[i], grid1=0; Temp; 
	 Temp = Temp->NextGridThisLevel, grid1++)
      if (UpdateMask[i][grid1]) {
      
	/* Get a list of possible siblings from the chaining mesh */
	
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

      } // ENDIF UpdateMask
    
#define NO_DEBUG
#ifdef DEBUG
    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel) {
      Temp->GridData->CheckSubgridMarker();
    }
#endif    

    FastSiblingLocatorFinalize(&ChainingMesh);      

  } // ENDFOR levels

  /* Cleanup */

  for (i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++) {
    if (LevelArray[i] != NULL)
      delete [] Grids[i];
    if (UpdateMask[i] != NULL)
      delete [] UpdateMask[i];
  }

  return SUCCESS;
}
