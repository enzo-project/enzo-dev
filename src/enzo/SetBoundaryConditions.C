/***********************************************************************
/
/  SET BOUNDARY CONDITIONS (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modifiedN:  Robert Harkness
/  date:       February, 2008
/
/  PURPOSE:
/======================================================================= 
/  This routine sets all the boundary conditions for Grids by either
/   interpolating from their parents or copying from sibling grids. 
/
/   This is part of a collection of routines called by EvolveLevel.
/   These have been optimized for enhanced message passing
/   performance by performing two passes -- one which generates
/   sends and the second which receives them.
/
/  modified: Robert Harkness, December 2007
/
************************************************************************/
 
#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"

/* function prototypes */
 
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

#define GRIDS_PER_LOOP 100000
 


#ifdef FAST_SIB
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry *Level)
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry *Level)
#endif
{
 
 
  int loopEnd = (ShearingBoundaryDirection != -1) ? 2 : 1;
  
 
  int grid1, grid2, StartGrid, EndGrid, loop;
  
  LCAPERF_START("SetBoundaryConditions");
  TIMER_START("SetBoundaryConditions");
    
  for (loop = 0; loop < loopEnd; loop++){
    
#ifdef FORCE_MSG_PROGRESS 
    CommunicationBarrier();
#endif
  
   

    if (loop == 0) {   

      TIME_MSG("Interpolating boundaries from parent");
  
    for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
	
      if (traceMPI) fprintf(tracePtr, "SBC loop\n");
	
      EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
	
      /* -------------- FIRST PASS ----------------- */
      /* Here, we just generate the calls to generate the receive buffers,
	 without actually doing anything. */
      
      LCAPERF_START("SetBC_Parent");
	
      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;
     
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
	
	/* a) Interpolate boundaries from the parent grid or set external
	   boundary conditions. */
	
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	

	  if (level == 0) {
    	    Grids[grid1]->GridData->SetExternalBoundaryValues(Exterior);
	  } else {
	    Grids[grid1]->GridData->InterpolateBoundaryFromParent
	      (Grids[grid1]->ParentGrid->GridData);
	  }
	

      } // ENDFOR grids

	/* -------------- SECOND PASS ----------------- */
	/* Now we generate all the sends, and do all the computation
	   for grids which are on the same processor as well. */

      CommunicationDirection = COMMUNICATION_SEND;
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

	/* a) Interpolate boundaries from the parent grid or set
	   external boundary conditions. */

	if (level > 0)
	  Grids[grid1]->GridData->InterpolateBoundaryFromParent
	    (Grids[grid1]->ParentGrid->GridData);
            
      }
      // ENDFOR grids

	//Grids[StartGrid]->GridData->PrintToScreenBoundaries(0);

	/* -------------- THIRD PASS ----------------- */
	/* In this final step, we get the messages as they come in and
	   then match them to the methods which generate the receive
	   handle. */

      if (CommunicationReceiveHandler() == FAIL)
	ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
      
    } // ENDFOR grid batches

    LCAPERF_STOP("SetBC_Parent");
    }
    TIME_MSG("Copying zones in SetBoundaryConditions");
    LCAPERF_START("SetBC_Siblings");
    for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
      EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

      /* -------------- FIRST PASS ----------------- */
      /* b) Copy any overlapping zones for sibling grids.  */

      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;
 
#ifdef FAST_SIB
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	  Grids[grid1]->GridData->
	    CheckForOverlap(SiblingList[grid1].GridList[grid2],
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyZonesFromGrid);
#else
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyZonesFromGrid);
#endif

      /* -------------- SECOND PASS ----------------- */
      /* b) Copy any overlapping zones for sibling grids.  */

      CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	  Grids[grid1]->GridData->
	    CheckForOverlap(SiblingList[grid1].GridList[grid2],
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyZonesFromGrid);
#else
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyZonesFromGrid);
#endif

      /* -------------- THIRD PASS ----------------- */

      if (CommunicationReceiveHandler() == FAIL)
	ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
      
    LCAPERF_STOP("SetBC_Siblings");

    } // end loop over batchs of grids

   // ENDFOR loop (for ShearingBox)
 
    /* c) Apply external reflecting boundary conditions, if needed.  */

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->CheckForExternalReflections
      (MetaData->LeftFaceBoundaryCondition,
       MetaData->RightFaceBoundaryCondition);

  
  
#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
  }
 
  TIMER_STOP("SetBoundaryConditions");
  LCAPERF_STOP("SetBoundaryConditions");

  return SUCCESS;
  
}
