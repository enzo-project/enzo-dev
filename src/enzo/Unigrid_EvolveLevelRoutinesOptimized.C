/***********************************************************************
/
/  EVOLVE LEVEL ROUTINES (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modified1:  Robert Harkness
/  date:       February, 2008
/
/  PURPOSE:  This is a collection of routines called by EvolveLevel.
/            These have been optimized for enhanced message passing
/            performance by performing two passes -- one which generates
/            sends and the second which receives them.
/
************************************************************************/

#ifdef USE_MPI 
#include <mpi.h>
#endif
#include <stdio.h>

#include "ErrorExceptions.h"
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
 
/* function prototypes */
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level,
                                FLOAT When);
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
                                   HierarchyEntry *Grids[], int NumberOfGrids);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
 
 
extern int CopyPotentialFieldAverage;
 
/* ======================================================================= */
/* This routine sets all the boundary conditions for Grids by either
   interpolating from their parents or copying from sibling grids. */


int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry *Level)
{
  

  int loopEnd=1;  if (ShearingBoundaryDirection != -1) loopEnd=2;
  
   for (int loop=0; loop < loopEnd; loop++){
  int grid, grid2;
 
  /* -------------- FIRST PASS ----------------- */

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND;
 
  if (traceMPI) fprintf(tracePtr, "SBC send\n");
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
 
    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */
 
    if (loop==0){
    if (level == 0) {
	if (Grids[grid]->GridData->SetExternalBoundaryValues(Exterior)
	  == FAIL) {
	ENZO_FAIL("Error in grid->SetExternalBoundaryValues.\n");
      
	}}
    else {
      if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid]->ParentGrid->GridData)) == FAIL) {
	ENZO_FAIL("Error in grid->InterpolateBoundaryFromParent.\n");
      }
    }}
 
    /* b) Copy any overlapping zones for sibling grids.  */
 
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	== FAIL) {
      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
      }

    /* c) Apply external reflecting boundary conditions, if needed.  */

    if (Grids[grid1]->GridData->CheckForExternalReflections(
							   MetaData->LeftFaceBoundaryCondition,
							   MetaData->RightFaceBoundaryCondition)
        == FAIL) {
      fprintf(stderr, "Error in grid->CheckForExternalReflections.\n");
    }
 
  } // end loop over grids

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  /* -------------- SECOND PASS ----------------- */
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
  if (traceMPI) fprintf(tracePtr, "SBC recv\n");
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
 
    /* a) Interpolate boundaries from the parent grid or set external
       boundary conditions. */
 
    if (loop==0){
    if (level > 0)
      if ((Grids[grid]->GridData->InterpolateBoundaryFromParent
	   (Grids[grid]->ParentGrid->GridData)) == FAIL) {
	ENZO_FAIL("Error in grid->InterpolateBoundaryFromParent.\n");
      }
    }

    /* b) Copy any overlapping zones for sibling grids.  */
 
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition,
				     &grid::CopyZonesFromGrid)
	== FAIL) {
      fprintf(stderr, "Error in grid->CopyZonesFromGrid.\n");
    }
 
  } // end loop over grids

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
  
  return SUCCESS;
  
 
   }
} 
 
 
/* ======================================================================= */
/* This routine prepares the density field for all the grids on this level,
   both particle and baryonic densities.  It also calculates the potential
   field if this is level 0 (since this involves communication). */

int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When)
{
 
  int grid, grid2;
 
  /* Set the time for evaluation of the fields, etc. */
 
  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
                   When*LevelArray[level]->GridData->ReturnTimeStep();
 
  /* If level is above MaximumGravityRefinementLevel, then just update the
     gravity at the MaximumGravityRefinementLevel. */
 
  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);
 
  /* Create an array (Grids) of all the grids. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
 
  /* Grids: Deposit particles in their GravitatingMassFieldParticles. */
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Send)\n");

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassField.\n");
    }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Receive)\n");
 
  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (DepositParticleMassField(Grids[grid], EvaluateTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassField.\n");
    }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  /* Grids: compute the GravitatingMassField (baryons & particles). */
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
				    level, When) == FAIL) {
      ENZO_FAIL("Error in PrepareGravitatingMassField.\n");
    }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (PrepareGravitatingMassField(Grids[grid], MetaData, LevelArray,
				    level, When) == FAIL) {
      ENZO_FAIL("Error in PrepareGravitatingMassField.\n");
    }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */
 
  //  if (level > 0)
 
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_SEND;

  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	ENZO_FAIL("Error in grid->CopyOverlappingMassField.\n");
      }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF2 (receive)\n", MyProcessorNumber);
 
  CommunicationDirection = COMMUNICATION_RECEIVE;

  for (grid = 0; grid < NumberOfGrids; grid++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyOverlappingMassField) == FAIL) {
	ENZO_FAIL("Error in grid->CopyOverlappingMassField.\n");
      }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /* Compute the potential for the top grid. */
 
  if (level == 0) {
  if (traceMPI) fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero (send-receive)\n", MyProcessorNumber);
    if (ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids) == FAIL) {
      ENZO_FAIL("Error in ComputePotentialFieldLevelZero.\n");
    }
  }
 
  /* Compute a first iteration of the potential and share BV's. */
 
#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
      if (level > 0) {
	CopyPotentialFieldAverage = 1;
	for (int iterate = 0; iterate < MAX_POTENTIAL_ITERATIONS; iterate++) {
 
	  if (iterate > 0)
	    CopyPotentialFieldAverage = 2;
 
	  int Dummy, grid2;
	  for (grid = 0; grid < NumberOfGrids; grid++) {
	    if (Grids[grid]->GridData->SolveForPotential(Dummy, level,
							 EvaluateTime)
		== FAIL) {
	      ENZO_FAIL("Error in grid->SolveForPotential.\n");
	    }
            if (CopyGravPotential)
            {
              // fprintf(stderr, "Call CP from PrepareDensityField 1\n");
              Grids[grid]->GridData->CopyPotentialToBaryonField();
            }
          }
 
          if (traceMPI) fprintf(tracePtr, "ITPOT send\n");

#ifdef FORCE_MSG_PROGRESS 
          CommunicationBarrier();
#endif

	  CommunicationDirection = COMMUNICATION_SEND;

	  for (grid = 0; grid < NumberOfGrids; grid++)
	    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	     if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       ENZO_FAIL("Error in grid->CopyPotentialField.\n");
	     }

#ifdef FORCE_MSG_PROGRESS 
          CommunicationBarrier();
#endif

          if (traceMPI) fprintf(tracePtr, "ITPOT recv\n");
 
	  CommunicationDirection = COMMUNICATION_RECEIVE;

	  for (grid = 0; grid < NumberOfGrids; grid++)
	    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	     if (Grids[grid]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				   MetaData->LeftFaceBoundaryCondition,
				   MetaData->RightFaceBoundaryCondition,
				   &grid::CopyPotentialField) == FAIL) {
	       ENZO_FAIL("Error in grid->CopyPotentialField.\n");
	     }

#ifdef FORCE_MSG_PROGRESS
          CommunicationBarrier();
#endif

	  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
	}
	CopyPotentialFieldAverage = 0;
      }
#endif /* ITERATE_POTENTIAL */
 
  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */
 
  if (reallevel > MaximumGravityRefinementLevel) {
 
    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
        if there is a subgrrid on reallevel, but this is ok). */
 
    for (grid = 0; grid < NumberOfGrids; grid++)
      if (Grids[grid]->NextGridNextLevel != NULL) {
	Grids[grid]->GridData->SolveForPotential(level,
					       MaximumGravityRefinementLevel);
        if (CopyGravPotential)
        {
          // fprintf(stderr, "Call CP from PrepareDensityField 2\n");
          Grids[grid]->GridData->CopyPotentialToBaryonField();
        }
        else
        {
          Grids[grid]->GridData->ComputeAccelerationField(
			(HydroMethod == Zeus_Hydro) ? DIFFERENCE_TYPE_STAGGERED : DIFFERENCE_TYPE_NORMAL,
                                               MaximumGravityRefinementLevel);
        }
      }
 
    /* Interpolate potential for reallevel grids from coarser grids. */
 
    if (!CopyGravPotential) {
 
      int Dummy;
      LevelHierarchyEntry *Temp = LevelArray[reallevel];

#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      CommunicationDirection = COMMUNICATION_SEND;

      while (Temp != NULL) {
        HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
        for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	  Temp3 = Temp3->ParentGrid;
        if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->InterpolateAccelerations.\n");
        }
        Temp = Temp->NextGridThisLevel;
      }

#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      CommunicationDirection = COMMUNICATION_RECEIVE;

      Temp = LevelArray[reallevel];
      while (Temp != NULL) {
        HierarchyEntry *Temp3 = Temp->GridHierarchyEntry;
        for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	  Temp3 = Temp3->ParentGrid;
        if (Temp->GridData->InterpolateAccelerations(Temp3->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->InterpolateAccelerations.\n");
        }
        Temp = Temp->NextGridThisLevel;
      }

#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
    } // end:  if (!CopyGravPotential)
 
  } // end: if (reallevel > MaximumGravityRefinementLevel)
 
  return SUCCESS;
}
 
 
/* ======================================================================= */
/* This routines does the flux correction and project for all grids on this
   level from the list of subgrids. */
 
int UpdateFromFinerGrids(int bogus_level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[])
 
{
 
  int grid, subgrid;
  HierarchyEntry *NextGrid;
 
  /* Define a temporary flux holder for the refined fluxes. */
 
  fluxes SubgridFluxesRefined;
 
  /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */
 
  /* -------------- FIRST PASS ----------------- */

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND;
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
 
    /* Loop over subgrids for this grid. */
 
    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {
 
      /* Project subgrid's refined fluxes to the level of this grid. */
 
      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	ENZO_FAIL("Error in grid->GetProjectedBoundaryFluxes.\n");
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }
 
    /* Loop over subgrids for this grid: replace solution. */
 
    NextGrid = Grids[grid]->NextGridNextLevel;
    while (NextGrid != NULL) {
 
      /* Project the subgrid solution into this grid. */
 
      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid]->GridData) == FAIL) {
	ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
    }
 
  } // end of loop over subgrids

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  /* -------------- SECOND PASS ----------------- */
 
  CommunicationDirection = COMMUNICATION_RECEIVE;
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
 
    /* Loop over subgrids for this grid. */
 
    NextGrid = Grids[grid]->NextGridNextLevel;
    subgrid = 0;
    while (NextGrid != NULL && FluxCorrection) {
 
      /* Project subgrid's refined fluxes to the level of this grid. */
 
      if (NextGrid->GridData->GetProjectedBoundaryFluxes(
		      Grids[grid]->GridData, SubgridFluxesRefined) == FAIL) {
	ENZO_FAIL("Error in grid->GetProjectedBoundaryFluxes.\n");
      }
	
      /* Correct this grid for the refined fluxes (step #19)
	 (this also deletes the fields in SubgridFluxesRefined). */
 
      if (Grids[grid]->GridData->CorrectForRefinedFluxes
	  (SubgridFluxesEstimate[grid][subgrid], &SubgridFluxesRefined,
	   SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1]     )
	  == FAIL) {
	ENZO_FAIL("Error in grid->CorrectForRefinedFluxes.\n");
      }
 
      NextGrid = NextGrid->NextGridThisLevel;
      subgrid++;
    }
 
    /* Loop over subgrids for this grid: replace solution. */
 
    NextGrid = Grids[grid]->NextGridNextLevel;
    while (NextGrid != NULL) {
 
      /* Project the subgrid solution into this grid. */
 
      if (NextGrid->GridData->ProjectSolutionToParentGrid
	                                   (*Grids[grid]->GridData) == FAIL) {
	ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");

      }
 
      NextGrid = NextGrid->NextGridThisLevel;
    }
 
  } // end of loop over subgrids

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
 
