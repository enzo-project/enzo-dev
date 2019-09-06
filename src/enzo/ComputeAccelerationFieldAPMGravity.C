/***********************************************************************
/
/  COMPUTE ACCELERATION FIELD FUNCTION
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

#include <stdio.h>
#include <stdlib.h>
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
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "communication.h"
#include "CommunicationUtilities.h"

/* function prototypes */

int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
                                int NumberOfSubgrids[] = NULL,
                                int FluxFlag = FALSE,
                                TopGridData* MetaData = NULL);

/* EvolveHierarchy function */

int ComputeAccelerationFieldAPMGravity(HierarchyEntry *Grid, TopGridData *MetaData,
			     LevelHierarchyEntry *LevelArray[], int level,
			     ExternalBoundary *Exterior = NULL)
{

  /* declarations */

  int dim;
  grid *CurrentGrid = Grid->GridData;

  /* Compute the refinement factor. */

  int RefinementFactor = 1, Factors[MAX_DIMENSION];
  if (Grid->ParentGrid != NULL) {
    Grid->ParentGrid->GridData->ComputeRefinementFactors(CurrentGrid, Factors);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (Factors[dim] != 1) 
	RefinementFactor = Factors[dim];
    }
  }

  /* Compute the acceleration (and potential) on the grid. */
  if (level == 0) {

    /* Here we calculate the acceleration field for PPM 
       which is assumed to be face-centered.
       If Zeus is used, staggering will be done in ZeusSource.C */
    int DiffType = DIFFERENCE_TYPE_NORMAL;
    if (CurrentGrid->ComputeAccelerationField(DiffType, level) == FAIL)
      ENZO_FAIL("Error in grid->ComputeAccelerationField.C\n");
    
  } else {
    if (CurrentGrid->ComputeAccelerationFieldAPM(RefinementFactor) == FAIL)
      ENZO_FAIL("Error in grid->ComputeAccelerationFieldAPM.\n");
  
    /* For the TestSelfForce we may want to output the force components separately */
    if (APMAddParentContribution) {

      
      /* a) Add parental acceleration. */
      /* Three-pass structure based on SetBoundaryConditions.C */
      LCAPERF_START("APMAddParentAcceleration");
      TIMER_START("APMAddParentAcceleration");

#ifdef FORCE_MSG_PROGRESS
      CommunicationBarrier();
#endif
      
      /* -------------- FIRST PASS ----------------- */
      /* Here, we just generate the calls to generate the receive buffers,
         without actually doing anything. */
      
      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;

      if (Grid->ParentGrid != NULL)
	if (CurrentGrid->AddParentAccelerationFieldAPM(Grid->ParentGrid->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->AddParentAccelerationFieldAPM.C\n");

      /* -------------- SECOND PASS ----------------- */
      /* Now we generate all the sends, and do all the computation
	 for grids which are on the same processor as well. */

      CommunicationDirection = COMMUNICATION_SEND;
      
      if (Grid->ParentGrid != NULL)
	if (CurrentGrid->AddParentAccelerationFieldAPM(Grid->ParentGrid->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->AddParentAccelerationFieldAPM.C\n");


      /* -------------- THIRD PASS ----------------- */
      /* In this final step, we get the messages as they come in and
	 then match them to the methods which generate the receive
	 handle. */

      if (CommunicationReceiveHandler() == FAIL)
        ENZO_FAIL("CommunicationReceiveHandler() failed!\n");


#ifdef FORCE_MSG_PROGRESS
      CommunicationBarrier();
#endif

      CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

      TIMER_STOP("APMAddParentAcceleration");
      LCAPERF_STOP("APMAddParentAcceleration");



      /* b) Add parental potential */
      /* Three-pass structure based on SetBoundaryConditions.C */
      LCAPERF_START("APMAddParentPotential");
      TIMER_START("APMAddParentPotential");

#ifdef FORCE_MSG_PROGRESS
      CommunicationBarrier();
#endif
      
      /* -------------- FIRST PASS ----------------- */
      /* Here, we just generate the calls to generate the receive buffers,
         without actually doing anything. */
      
      CommunicationDirection = COMMUNICATION_POST_RECEIVE;
      CommunicationReceiveIndex = 0;

      if (Grid->ParentGrid != NULL)
	if (CurrentGrid->AddParentPotentialFieldAPM(Grid->ParentGrid->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->AddParentAccelerationFieldAPM.C\n");

      /* -------------- SECOND PASS ----------------- */
      /* Now we generate all the sends, and do all the computation
	 for grids which are on the same processor as well. */

      CommunicationDirection = COMMUNICATION_SEND;
      
      if (Grid->ParentGrid != NULL)
	if (CurrentGrid->AddParentPotentialFieldAPM(Grid->ParentGrid->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->AddParentAccelerationFieldAPM.C\n");


      /* -------------- THIRD PASS ----------------- */
      /* In this final step, we get the messages as they come in and
	 then match them to the methods which generate the receive
	 handle. */

      if (CommunicationReceiveHandler() == FAIL)
        ENZO_FAIL("CommunicationReceiveHandler() failed!\n");


#ifdef FORCE_MSG_PROGRESS
      CommunicationBarrier();
#endif

      CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

      TIMER_STOP("APMAddParentPotential");
      LCAPERF_STOP("APMAddParentPotential");

    } // end if APMAddParentContribution
  }/* end if (level == 0) */
  
  /* Clear grid and particle accelerations. */
  
  //  CurrentGrid->ComputeEffectiveGridPositions(0.0, TRUE, 1.0);
  //  CurrentGrid->ClearAccelerationFieldForCells();
  CurrentGrid->ClearParticleAccelerations();

  /* Move particles 1/2 step forward in preparation for interpolation. */

  CurrentGrid->UpdateParticlePosition(0.5*CurrentGrid->ReturnTimeStep());

  /* Zeus: reset grid positions to cell centers (move to faces later). */

  //    if (HydroMethod == Zeus_Hydro)
  //      CurrentGrid->ComputeEffectiveGridPositions(0.0, TRUE, 1.0);

  /* Interpolate the accelerations back to the grid and particles. */

  if (level <= MaximumGravityRefinementLevel) {
    CurrentGrid->InterpolateParticlePositions(CurrentGrid, DIFFERENCE_TYPE_NORMAL);

    //    if (Grid->ParentGrid != NULL)
    //      CurrentGrid->InterpolateParticlePositions(Grid->ParentGrid->GridData, DIFFERENCE_TYPE_NORMAL);
    //    else
    //      CurrentGrid->InterpolateParticlePositions(CurrentGrid, DIFFERENCE_TYPE_NORMAL);
  }

  /* Interpolate the accelerations of all grids higher than this one
     to the particles and grids. */

#ifdef UNUSED
  HierarchyEntry *Temp = Grid->ParentGrid;
  int l1 = level-1;
  while (Temp != NULL) {
    if (l1-- <= MaximumGravityRefinementLevel) {
      //	if (CurrentGrid->InterpolateGridPositions(Temp->GridData) == FAIL) {
      //	  fprintf(stderr, "Error in grid->InterpolateGridPositions.\n");
      //	  return FAIL;
      //	}
      if (CurrentGrid->InterpolateParticlePositions(Temp->GridData, DIFFERENCE_TYPE_NORMAL) == FAIL)
	ENZO_FAIL("Error in grid->InterpolateParticlePositons.\n");
    }
    Temp = Temp->ParentGrid;
  }
#endif

  /* Copy the cell accelerations to a grid.  This routine really just
     readjusts the size of AccelerationFieldForCells so that it is
     commensurate with the rest of the baryonic grids. */

#if 0
  if (CurrentGrid->CopyGridAccelerationsToGrid() == FAIL)
    ENZO_FAIL("Error in grid->CopyGridAccelerationsToGrid.\n");
#endif

  /* Move particles 1/2 step backwards to return to original positions. */

  CurrentGrid->UpdateParticlePosition(-0.5*CurrentGrid->ReturnTimeStep());

  /* Set any boundary conditions needed (only for top grid).
     For now, we'll use the same boundary conditions as the density.  This is
     correct for periodic B.C.'s but probably not for open B.C.'s. */
  // FIX

#if 0
  if (Grid->ParentGrid == NULL)
    if (CurrentGrid->SetGravitationalExternalBoundary(Exterior) == FAIL)
      ENZO_FAIL("Error in grid->SetGravitationalExternalBoundary.\n");
#endif 

  //  CurrentGrid->DeleteGravitatingMassField();
  //  CurrentGrid->DeleteGridPositionAndMass();
  
  /* If this is the lowest level of the hierarchy, we can delete
     AccelerationField as well. */

  //  if (LevelArray[level+1] == NULL)
  //    CurrentGrid->DeleteAccelerationField();

  return SUCCESS;
}
