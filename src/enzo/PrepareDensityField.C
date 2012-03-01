/***********************************************************************
/
/  PREPARE DENSITY FIELD (CALLED BY EVOLVE LEVEL)
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modifiedN:  Robert Harkness
/  date:       February, 2008
/
/ ======================================================================= 
/ This routine prepares the density field for all the grids on this level,
/ both particle and baryonic densities.  It also calculates the potential
/ field if this is level 0 (since this involves communication). 
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
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
 
int PrepareGravitatingMassField1(HierarchyEntry *Grid);
#ifdef FAST_SIB
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When);
#else
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When);
#endif

int PrepareGravitatingMassField2b(HierarchyEntry *Grid, int level);
 
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
 
 
 
extern int CopyPotentialFieldAverage;
 
#define GRIDS_PER_LOOP 100000

 
#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When)
#else   // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When)
#endif  // end FAST_SIB
{

  /* Return if this does not concern us */
  if (!SelfGravity) return SUCCESS;
 
  LCAPERF_START("PrepareDensityField");

  int grid1, grid2, StartGrid, EndGrid;
 
  /* Set the time for evaluation of the fields, etc. */
 
  FLOAT EvaluateTime = LevelArray[level]->GridData->ReturnTime() +
    When*LevelArray[level]->GridData->ReturnTimeStep();
 
  /* If level is above MaximumGravityRefinementLevel, then just
     update the gravity at the MaximumGravityRefinementLevel. */
 
  int reallevel = level;
  level = min(level, MaximumGravityRefinementLevel);
 
  /* Create an array (Grids) of all the grids. */
 
  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  /************************************************************************/
  /* Grids: Deposit particles in their GravitatingMassFieldParticles.
     (Do a batch of grids at a time; this is a loop over the batches)
  */

  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField (Send)\n");

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  TIME_MSG("Depositing particle mass field");
  LCAPERF_START("DepositParticleMassField");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      DepositParticleMassField(Grids[grid1], EvaluateTime);

#ifdef FORCE_MSG_PROGRESS 
    CommunicationBarrier();
#endif

    if (traceMPI) 
      fprintf(tracePtr, "PrepareDensityField: Enter DepositParticleMassField"
	      " (Receive)\n");
 
    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      DepositParticleMassField(Grids[grid1], EvaluateTime);

    /* Finally, receive the data and process it. */
    
    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("DepositParticleMassField");
    

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  /******************************************************************/
  /* Grids: compute the GravitatingMassField (baryons & particles). */
  /*   This is now split into two section. */
 
  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF1 (send)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("PrepareGravitatingMassField1");
  LCAPERF_START("PrepareGravitatingMassField1");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 1 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
 
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField1(Grids[grid1]);

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField1(Grids[grid1]);

    /* Finally, receive the data and process it. */
    
    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField1");


#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("PrepareGravitatingMassField2");
  LCAPERF_START("PrepareGravitatingMassField2a");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef BITWISE_IDENTICALITY
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
      
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
				    MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
				    level, When);
#endif

#ifndef BITWISE_IDENTICALITY
    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], grid1, SiblingList,
				   MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2a(Grids[grid1], MetaData, LevelArray,
				   level, When);
#endif

    CommunicationReceiveHandler();
#endif /* BITWISE_IDENTICALITY */

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField2a");

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  /************************************************************************/
  LCAPERF_START("PrepareGravitatingMassField2b");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2b(Grids[grid1], level);

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2b(Grids[grid1], level);

    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("PrepareGravitatingMassField2b");

  /************************************************************************/
  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */
 
  //  if (level > 0)
 
  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("CopyOverlappingMassField");
  LCAPERF_START("CopyOverlappingMassField");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
      
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(SiblingList[grid1].GridList[grid2],

			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(Grids[grid2]->GridData,
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#endif

    CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(SiblingList[grid1].GridList[grid2],
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	Grids[grid1]->GridData->
	  CheckForOverlap(Grids[grid2]->GridData,
			  MetaData->LeftFaceBoundaryCondition,
			  MetaData->RightFaceBoundaryCondition,
			  &grid::CopyOverlappingMassField);
#endif

    CommunicationReceiveHandler();

  } // ENDFOR grid batches
  LCAPERF_STOP("CopyOverlappingMassField");

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /************************************************************************/
  /* Compute the potential for the top grid. */
 
  if (level == 0) {
    TIME_MSG("ComputePotentialFieldLevelZero");
    LCAPERF_START("ComputePotentialFieldLevelZero");
    TIMER_START("ComputePotentialFieldLevelZero");
    if (traceMPI) 
      fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero "
	      "(send-receive)\n", MyProcessorNumber);
#ifdef FAST_SIB
    ComputePotentialFieldLevelZero(MetaData, SiblingList,
				   Grids, NumberOfGrids);
#else
    ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids);
#endif
    TIMER_STOP("ComputePotentialFieldLevelZero");
    LCAPERF_STOP("ComputePotentialFieldLevelZero");
  }
       
  /************************************************************************/
  /* Compute a first iteration of the potential and share BV's. */
 
  int iterate;
  if (level > 0) {
    LCAPERF_START("SolveForPotential");
    TIMER_START("SolveForPotential");
    CopyPotentialFieldAverage = 1;
    for (iterate = 0; iterate < PotentialIterations; iterate++) {
      
      if (iterate > 0)
	CopyPotentialFieldAverage = 2;

 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
	Grids[grid1]->GridData->SolveForPotential(level, EvaluateTime);
	if (CopyGravPotential)
	  Grids[grid1]->GridData->CopyPotentialToBaryonField();
      }
 
      if (traceMPI) fprintf(tracePtr, "ITPOT post-recv\n");
	
#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      TIME_MSG("CopyPotentialField");
      for (StartGrid = 0; StartGrid < NumberOfGrids; 
	   StartGrid += GRIDS_PER_LOOP) {
	EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);
  
#ifdef BITWISE_IDENTICALITY
	CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
#else
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
#endif
	CommunicationReceiveIndex = 0;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
#ifdef FAST_SIB
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
 
	  //fprintf(stderr, "#SIBSend on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
	  // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(SiblingList[grid1].GridList[grid2],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
	    
	  grid2 = grid1;
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyPotentialField);
	  
	} // ENDFOR grid1
#else
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	  for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(Grids[grid2]->GridData,
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
#endif

#ifndef BITWISE_IDENTICALITY
#ifdef FORCE_MSG_PROGRESS 
	CommunicationBarrier();
#endif

	if (traceMPI) fprintf(tracePtr, "ITPOT send\n");
 
	CommunicationDirection = COMMUNICATION_SEND;

 
#ifdef FAST_SIB
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++) {
 
	  //fprintf(stderr, "#SIBRecv on cpu %"ISYM": %"ISYM"\n", MyProcessorNumber, SiblingList[grid1].NumberOfSiblings);
 
	  // for (grid2 = SiblingList[grid1].NumberOfSiblings-1; grid2 = 0; grid2--)
	  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(SiblingList[grid1].GridList[grid2],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
 
	  grid2 = grid1;
	  Grids[grid1]->GridData->
	    CheckForOverlap(Grids[grid2]->GridData,
			    MetaData->LeftFaceBoundaryCondition,
			    MetaData->RightFaceBoundaryCondition,
			    &grid::CopyPotentialField);
 
	} // ENDFOR grid1
#else
	for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
	  for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
	    Grids[grid1]->GridData->
	      CheckForOverlap(Grids[grid2]->GridData,
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition,
			      &grid::CopyPotentialField);
#endif

	CommunicationReceiveHandler();
#endif

      } // ENDFOR grid batches
    } // ENDFOR iterations
    CopyPotentialFieldAverage = 0;
    TIMER_STOP("SolveForPotential");
    LCAPERF_STOP("SolveForPotential");
  } // ENDIF level > 0
  
  /* if level > MaximumGravityRefinementLevel, then do final potential
     solve (and acceleration interpolation) here rather than in the main
     EvolveLevel since it involves communications. */
  
  if (reallevel > MaximumGravityRefinementLevel) {
 
    /* compute potential and acceleration on coarser level [LOCAL]
       (but only if there is at least a subgrid -- it should be only
       if there is a subgrrid on reallevel, but this is ok). */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      if (Grids[grid1]->NextGridNextLevel != NULL) {
	Grids[grid1]->GridData->SolveForPotential(MaximumGravityRefinementLevel);
	if (CopyGravPotential)
	  Grids[grid1]->GridData->CopyPotentialToBaryonField();
	else
	  Grids[grid1]->GridData->ComputeAccelerationField
	    ((HydroMethod == Zeus_Hydro) ? DIFFERENCE_TYPE_STAGGERED : 
	     DIFFERENCE_TYPE_NORMAL, MaximumGravityRefinementLevel);
      }
 
    /* Interpolate potential for reallevel grids from coarser grids. */
 
    if (!CopyGravPotential) {
 
      int Dummy, GridCount;
      LevelHierarchyEntry *Temp, *LastTemp;
      HierarchyEntry *Temp3;
      LevelHierarchyEntry *FirstTemp = LevelArray[reallevel];
	
#ifdef FORCE_MSG_PROGRESS 
      CommunicationBarrier();
#endif

      do {

	GridCount = 0;
	CommunicationDirection = COMMUNICATION_POST_RECEIVE;
	CommunicationReceiveIndex = 0;
	CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
	Temp = FirstTemp;
	while (Temp != NULL && GridCount++ < GRIDS_PER_LOOP) {
	  Temp3 = Temp->GridHierarchyEntry;
	  for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	    Temp3 = Temp3->ParentGrid;
	  Temp->GridData->InterpolateAccelerations(Temp3->GridData);
	  Temp = Temp->NextGridThisLevel;
	} // ENDWHILE
	LastTemp = Temp;

	CommunicationDirection = COMMUNICATION_SEND;
	Temp = FirstTemp;
	while (Temp != LastTemp) {
	  Temp3 = Temp->GridHierarchyEntry;
	  for (Dummy = reallevel; Dummy > MaximumGravityRefinementLevel; Dummy--)
	    Temp3 = Temp3->ParentGrid;
	  Temp->GridData->InterpolateAccelerations(Temp3->GridData);
	  Temp = Temp->NextGridThisLevel;
	}
	FirstTemp = LastTemp;

	CommunicationReceiveHandler();

      } while (LastTemp != NULL);

    } // end:  if (!CopyGravPotential)
 
  } // end: if (reallevel > MaximumGravityRefinementLevel)

    // --------------------------------------------------
    // MEMORY LEAK FIX
    //
    // valgrind error: "1,388,304 (67,352 direct, 1,320,952 indirect)
    // bytes in 130 blocks are definitely lost in loss record 22 of 46"
    //
    // Adding missing delete [] () for Grids[] allocated in
    // GenerateGridArray()
    // --------------------------------------------------

  delete [] Grids;

  // --------------------------------------------------

  LCAPERF_STOP("PrepareDensityField");
  return SUCCESS;

}
 
