/***********************************************************************
/
/  COPY ZONES (non-blocking) FROM OLD TO NEW GRIDS (REBUILD HIERARCHY)
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
 
#include <stdio.h>
#include <string.h>
#include <time.h>
 
#include "ErrorExceptions.h"
#include "performance.h"
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
#include "CommunicationUtilities.h"
#include "communication.h"

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int CopyZonesFromOldGrids(LevelHierarchyEntry *OldGrids, 
			  TopGridData *MetaData,
			  ChainingMeshStructure ChainingMesh)
{

  int i;
  FLOAT ZeroVector[] = {0,0,0};
  LevelHierarchyEntry *Temp;
  SiblingGridList SiblingList;

  /* Post receives, looping over the old grids */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel) {

    // Find sibling grids
    Temp->GridData->FastSiblingLocatorFindSiblings
      (&ChainingMesh, &SiblingList, MetaData->LeftFaceBoundaryCondition, 
       MetaData->RightFaceBoundaryCondition);

    // For each of the sibling grids, copy data.
    for (i = 0; i < SiblingList.NumberOfSiblings; i++)
      SiblingList.GridList[i]->CopyZonesFromGrid(Temp->GridData, ZeroVector);

    // Don't delete the old grids yet, we need to copy their data in
    // the next step.
    
    delete [] SiblingList.GridList;

  }

  /* Send data */

  CommunicationDirection = COMMUNICATION_SEND;

  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel) {

    /* Find sibling grids.  Note that we do this twice to save memory.
       This step isn't that expensive, so we can afford to compute
       this twice.  However this may change in larger simulations. */

    Temp->GridData->FastSiblingLocatorFindSiblings
      (&ChainingMesh, &SiblingList, MetaData->LeftFaceBoundaryCondition, 
       MetaData->RightFaceBoundaryCondition);

    // For each of the sibling grids, copy data.
    for (i = 0; i < SiblingList.NumberOfSiblings; i++)
      SiblingList.GridList[i]->CopyZonesFromGrid(Temp->GridData, ZeroVector);

    /* Delete all fields (only on the host processor -- we need
       BaryonField on the receiving processor) after sending them.  We
       only delete the grid object on all processors after
       everything's done. */

    if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      Temp->GridData->DeleteAllFields();

    delete [] SiblingList.GridList;

  }

  /* Receive data */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("");

  /* Delete old grids */

  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel) {
    delete Temp->GridData;
    Temp->GridData = NULL;
  }

  return SUCCESS;

}
