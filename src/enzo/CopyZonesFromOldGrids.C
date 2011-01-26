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

#define GRIDS_PER_LOOP 100000
#define CELLS_PER_LOOP 100000000

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int CopyZonesFromOldGrids(LevelHierarchyEntry *OldGrids, 
			  TopGridData *MetaData,
			  ChainingMeshStructure ChainingMesh)
{

  int i, dim, size, NumberOfGrids, gridcount, totalcount, Rank, ncells;
  int EndGrid, Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  FLOAT ZeroVector[] = {0,0,0};
  LevelHierarchyEntry *Temp, *FirstGrid;
  SiblingGridList SiblingList;

  /* Count grids */

  NumberOfGrids = 0;
  for (Temp = OldGrids; Temp; Temp = Temp->NextGridThisLevel)
    NumberOfGrids++;

  /* Loop over batches of grids */

  totalcount = 0;
  FirstGrid = OldGrids;
  while (totalcount < NumberOfGrids && FirstGrid != NULL) {

    /* Find how many grids have CELLS_PER_LOOP */

    for (Temp = FirstGrid, ncells = 0, gridcount = 0;
	 Temp && ncells < CELLS_PER_LOOP && gridcount < GRIDS_PER_LOOP;
	 Temp = Temp->NextGridThisLevel, gridcount++) {
      Temp->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
      size = 1;
      for (dim = 0; dim < Rank; dim++)
	size *= Dims[dim];
      ncells += size;
    }

    EndGrid = gridcount;

    /* Post receives, looping over the old grids */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    
    for (Temp = FirstGrid, gridcount = 0; 
	 Temp && gridcount < EndGrid; 
	 Temp = Temp->NextGridThisLevel, gridcount++) {

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

    for (Temp = FirstGrid, gridcount = 0; 
	 Temp && gridcount < EndGrid; 
	 Temp = Temp->NextGridThisLevel, gridcount++) {

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
      ENZO_FAIL("CommunicationReceiveHandler() failed!\n");


    /* Delete old grids and increase total grid count and then advance
       FirstGrid pointer */

    for (Temp = FirstGrid, gridcount = 0; 
	 Temp && gridcount < EndGrid; 
	 Temp = Temp->NextGridThisLevel, gridcount++) {
      delete Temp->GridData;
      Temp->GridData = NULL;
      totalcount++;
    }

    FirstGrid = Temp;
#ifdef USE_MPI
    CommunicationBufferPurge();
#endif /* USE_MPI */
  } // ENDWHILE grid batches

  return SUCCESS;

}
