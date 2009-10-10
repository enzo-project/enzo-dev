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
 
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

#define GRIDS_PER_LOOP 20000
 
 
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
 
  
