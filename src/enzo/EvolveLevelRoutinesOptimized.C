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
 
 
 
extern int CopyPotentialFieldAverage;
 
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
  
  for (loop = 0; loop < loopEnd; loop++){
    
    JBPERF_START("SetBoundaryConditions");
    
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

  } // ENDFOR loop
 
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
  /* This routine prepares the density field for all the grids on this level,
     both particle and baryonic densities.  It also calculates the potential
     field if this is level 0 (since this involves communication). */
 
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
 
  JBPERF_START("PrepareDensityField");

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


#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): PGMF2 (receive)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("PrepareGravitatingMassField2");
  for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    EndGrid = min(StartGrid + GRIDS_PER_LOOP, NumberOfGrids);

    /* ----- section 2 ---- */
    /* First, generate the receive calls. */

    CommunicationReceiveIndex = 0;
    CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
    CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2(Grids[grid1], grid1, SiblingList,
				   MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2(Grids[grid1], MetaData, LevelArray,
				   level, When);
#endif

    /* Next, send data and process grids on the same processor. */

    CommunicationDirection = COMMUNICATION_SEND;
#ifdef FAST_SIB
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2(Grids[grid1], grid1, SiblingList,
				   MetaData, level, When);
#else
    for (grid1 = StartGrid; grid1 < EndGrid; grid1++)
      PrepareGravitatingMassField2(Grids[grid1], MetaData, LevelArray,
				   level, When);
#endif

    CommunicationReceiveHandler();

  } // ENDFOR grid batches

#ifdef FORCE_BUFFER_PURGE
  CommunicationBufferPurge();
#endif

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif
 
  /************************************************************************/
  /* Copy overlapping mass fields to ensure consistency and B.C.'s. */
 
  //  if (level > 0)
 
  if (traceMPI) 
    fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): COMF1 (send)\n", 
	    MyProcessorNumber);
 
  TIME_MSG("CopyOverlappingMassField");
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
    if (traceMPI) 
      fprintf(tracePtr, "PrepareDensityField: P(%"ISYM"): CPFLZero "
	      "(send-receive)\n", MyProcessorNumber);
#ifdef FAST_SIB
    ComputePotentialFieldLevelZero(MetaData, SiblingList,
				   Grids, NumberOfGrids);
#else
    ComputePotentialFieldLevelZero(MetaData, Grids, NumberOfGrids);
#endif
  }
       
  /************************************************************************/
  /* Compute a first iteration of the potential and share BV's. */
 
#define ITERATE_POTENTIAL
#ifdef ITERATE_POTENTIAL
  int iterate;
  if (level > 0) {
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
  
	CommunicationDirection = COMMUNICATION_POST_RECEIVE;
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

      } // ENDFOR grid batches
    } // ENDFOR iterations
    CopyPotentialFieldAverage = 0;
  } // ENDIF level > 0
#endif /* ITERATE_POTENTIAL */
  
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

  JBPERF_STOP("PrepareDensityField");
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
 
  
