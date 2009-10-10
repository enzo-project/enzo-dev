/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modifiedN:  Robert Harkness
/  date:       February, 2008
/
/  PURPOSE:  This is a collection of routines called by EvolveLevel.
/            These have been optimized for enhanced message passing
/            performance by performing two passes -- one which generates
/            sends and the second which receives them.
/
/  modified: Robert Harkness, December 2007
/
************************************************************************/
 
#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */
 
#include <stdio.h>
#include "ErrorExceptions.h"
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
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
 
int PrepareGravitatingMassField1(HierarchyEntry *Grid);
#ifdef FAST_SIB
int PrepareGravitatingMassField2(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When);
#else
int PrepareGravitatingMassField2(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When);
#endif
 
#ifdef FAST_SIB
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   SiblingGridList SiblingList[],
				   HierarchyEntry *Grids[], int NumberOfGrids);
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], int NumberOfGrids);
#endif

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
 
 

#define GRIDS_PER_LOOP 20000
 
/* ======================================================================= */
/* This routine sets all the boundary conditions for Grids by either
   interpolating from their parents or copying from sibling grids. */
 


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
 
  int loopEnd = (ShearingBoundaryDirection != 1) ? 2 : 1;
  
  int grid1, grid2, StartGrid, EndGrid, loop;
  
  JBPERF_START("SetBoundaryConditions");
    
  for (loop = 0; loop < loopEnd; loop++){
    
#ifdef FORCE_MSG_PROGRESS 
    CommunicationBarrier();
#endif
  
    TIME_MSG("Interpolating boundaries from parent");
      
    for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
	
      if (traceMPI) fprintf(tracePtr, "SBC loop\n");
	
      EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
	
      /* -------------- FIRST PASS ----------------- */
      /* Here, we just generate the calls to generate the receive buffers,
	 without actually doing anything. */
      
      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;
      
      for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
	
	/* a) Interpolate boundaries from the parent grid or set external
	   boundary conditions. */
	
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	
	if (level == 0) {
	  if (loop == 0) 
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
     
      
      } // ENDFOR grids

	//Grids[StartGrid]->GridData->PrintToScreenBoundaries(0);

	/* -------------- THIRD PASS ----------------- */
	/* In this final step, we get the messages as they come in and
	   then match them to the methods which generate the receive
	   handle. */

      if (CommunicationReceiveHandler() == FAIL)
	ENZO_FAIL("");

    } // ENDFOR grid batches

    TIME_MSG("Copying zones in SetBoundaryConditions");
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
	ENZO_FAIL("");

    } // end loop over batchs of grids

  } // ENDFOR loop (for ShearinBox)
 
    /* c) Apply external reflecting boundary conditions, if needed.  */

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->CheckForExternalReflections
      (MetaData->LeftFaceBoundaryCondition,
       MetaData->RightFaceBoundaryCondition);
  
#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  JBPERF_STOP("SetBoundaryConditions");

  return SUCCESS;
  
}
 
  /* ======================================================================= */
  /* This routines does the flux correction and project for all grids on this
     level from the list of subgrids. */
 
#ifdef FLUX_FIX
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry* SUBlingList[],
			 TopGridData *MetaData)
#else
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[])
#endif
 
{
 
  int grid1, subgrid, StartGrid, EndGrid;
  HierarchyEntry *NextGrid;
 
#ifdef FLUX_FIX
  int SUBlingGrid;
  LevelHierarchyEntry *NextEntry;
#endif
 
  /* Define a temporary flux holder for the refined fluxes. */
 
  fluxes SubgridFluxesRefined;
 
  /* For each grid,
     (a)  project the subgrid's solution into this grid (step #18)
     (a1) project the neighboring subgrid's solution into this grid
     (b)  correct for the difference between this grid's fluxes and the
     subgrid's fluxes. (step #19) */
  
#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  TIME_MSG("UpdateFromFinerGrids");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* -------------- FIRST PASS ----------------- */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
      
      /* Loop over subgrids for this grid. */
 
      NextGrid = Grids[grid1]->NextGridNextLevel;
      subgrid = 0;
      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      
      while (NextGrid != NULL && FluxCorrection) {
 
	/* Project subgrid's refined fluxes to the level of this grid. */
 
#ifdef USE_MPI
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = grid1;
	CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = subgrid;
#endif /* USE_MPI */

	NextGrid->GridData->
	  GetProjectedBoundaryFluxes(Grids[grid1]->GridData, SubgridFluxesRefined);
 
	NextGrid = NextGrid->NextGridThisLevel;
	subgrid++;
      } // ENDWHILE subgrids

    } // ENDFOR grids

    /* -------------- SECOND PASS ----------------- */

    CommunicationDirection = COMMUNICATION_SEND;

    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid. */

      NextGrid = Grids[grid1]->NextGridNextLevel;
      subgrid = 0;
      while (NextGrid != NULL && FluxCorrection) {

	/* Project subgrid's refined fluxes to the level of this grid. */

	NextGrid->GridData->
	  GetProjectedBoundaryFluxes(Grids[grid1]->GridData, SubgridFluxesRefined);

	/* Correct this grid for the refined fluxes (step #19)
	   (this also deletes the fields in SubgridFluxesRefined). 
	   (only call it if the grid and sub-grid are on the same
	   processor, otherwise handled in CommunicationReceiveHandler.) */

#ifdef FLUX_FIX	
	if (NextGrid->GridData->ReturnProcessorNumber() ==
	    Grids[grid1]->GridData->ReturnProcessorNumber())
	  Grids[grid1]->GridData->CorrectForRefinedFluxes
	    (SubgridFluxesEstimate[grid1][subgrid], &SubgridFluxesRefined,
	     SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
	     FALSE, MetaData);
#else
	if (NextGrid->GridData->ReturnProcessorNumber() ==
	    Grids[grid1]->GridData->ReturnProcessorNumber())
	  Grids[grid1]->GridData->CorrectForRefinedFluxes
	    (SubgridFluxesEstimate[grid1][subgrid], &SubgridFluxesRefined,
	     SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1]);
#endif /* FLUX_FIX */

	NextGrid = NextGrid->NextGridThisLevel;
	subgrid++;

      } // ENDWHILE subgrids

    } // ENDFOR grids

    /* -------------- THIRD PASS ----------------- */

#ifdef FLUX_FIX
    CommunicationReceiveHandler(SubgridFluxesEstimate, NumberOfSubgrids, 
				FALSE, MetaData);
#else
    CommunicationReceiveHandler(SubgridFluxesEstimate, NumberOfSubgrids);
#endif

  } // ENDFOR grid batches

  /************************************************************************
     (a1) project the neighboring subgrid's solution into this grid 
          (SUBlings)
  ************************************************************************/

#ifdef FLUX_FIX
  TIME_MSG("Projecting neighboring subgrid solution (FLUX_FIX)");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* -------------- FIRST PASS ----------------- */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
	  
      /* Loop over subgrids for this grid. */
 
      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      NextEntry = SUBlingList[grid1];

      while (NextEntry != NULL && FluxCorrection) {

	/* make sure this isn't a "proper" subgrid */
	
	if (NextEntry->GridHierarchyEntry->ParentGrid != Grids[grid1]) {
	  
	  /* Project subgrid's refined fluxes to the level of this grid. */

#ifdef USE_MPI
	  CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = grid1;
	  CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = 
	    NumberOfSubgrids[grid1]-1;
#endif /* USE_MPI */

	  NextEntry->GridData->GetProjectedBoundaryFluxes
	    (Grids[grid1]->GridData, SubgridFluxesRefined);
	}

	NextEntry = NextEntry->NextGridThisLevel;

      } // ENDWHILE subgrids
    } // ENDFOR grids

    /* -------------- SECOND PASS ----------------- */

    CommunicationDirection = COMMUNICATION_SEND;

    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      NextEntry = SUBlingList[grid1];
 
      while (NextEntry != NULL && FluxCorrection) {

	/* make sure this isn't a "proper" subgrid */

	if (NextEntry->GridHierarchyEntry->ParentGrid != Grids[grid1]) {
	  
	  /* Project subgrid's refined fluxes to the level of this grid. */

	  NextEntry->GridData->GetProjectedBoundaryFluxes
	    (Grids[grid1]->GridData, SubgridFluxesRefined);
 
	  /* Correct this grid for the refined fluxes (step #19)
	     (this also deletes the fields in SubgridFluxesRefined). */
 
	  if (NextEntry->GridData->ReturnProcessorNumber() ==
	      Grids[grid1]->GridData->ReturnProcessorNumber())
	    Grids[grid1]->GridData->CorrectForRefinedFluxes
	      (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
	       &SubgridFluxesRefined,
	       SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1],
	       TRUE, MetaData);
	}

	NextEntry = NextEntry->NextGridThisLevel;

      } // ENDWHILE subgrids
    } // ENDFOR grids

    /* -------------- THIRD PASS ----------------- */

    CommunicationReceiveHandler(SubgridFluxesEstimate, NumberOfSubgrids, 
				TRUE, MetaData);

  } // ENDFOR grid batches
#endif /* FLUX_FIX */

  /************************************************************************
    (b) correct for the difference between this grid's fluxes and the
        subgrid's fluxes. (step #19) 
  ************************************************************************/

  TIME_MSG("Projecting solution to parent");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* -------------- FIRST PASS ----------------- */

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid: replace solution. */

      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      NextGrid = Grids[grid1]->NextGridNextLevel;
      while (NextGrid != NULL) {

	/* Project the subgrid solution into this grid. */

	NextGrid->GridData->ProjectSolutionToParentGrid(*Grids[grid1]->GridData);
	NextGrid = NextGrid->NextGridThisLevel;
      } // ENDWHILE subgrids
    } // ENDFOR grids

    /* -------------- SECOND PASS ----------------- */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {

      /* Loop over subgrids for this grid: replace solution. */

      CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      NextGrid = Grids[grid1]->NextGridNextLevel;
      while (NextGrid != NULL) {

	/* Project the subgrid solution into this grid. */

	NextGrid->GridData->ProjectSolutionToParentGrid(*Grids[grid1]->GridData);
	NextGrid = NextGrid->NextGridThisLevel;
      } // ENDWHILE subgrids
    } // ENDFOR grids

    /* -------------- THIRD PASS ----------------- */

    CommunicationReceiveHandler();

  } // ENDFOR grid batches


#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  return SUCCESS;
}
 
 
 
/* ======================================================================= */
/* This routine simply converts a linked list of grids into an array of
   pointers. */
 
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[])
{
 
  /* Count the number of grids on this level. */
 
  int NumberOfGrids = 0, counter = 0;
  LevelHierarchyEntry *Temp = LevelArray[level];
  while (Temp != NULL) {
    NumberOfGrids++;
    Temp             = Temp->NextGridThisLevel;
  }
 
  /* Create a list of pointers and number of subgrids (and fill it out). */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  *Grids = new HierarchyEntryPointer[NumberOfGrids];
  Temp = LevelArray[level];
  while (Temp != NULL) {
    (*Grids)[counter++] = Temp->GridHierarchyEntry;
    Temp              = Temp->NextGridThisLevel;
  }
 
  return NumberOfGrids;
}
 
  
