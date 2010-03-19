/***********************************************************************
/
/  PREPARE THE EMISSIVITY FIELD AND CALL FLD SOLVER
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/ PURPOSE: 
/          
/
************************************************************************/
#include "preincludes.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "StarParticleData.h"
#include "RadiativeTransferHealpixRoutines.h"
#include "ImplicitProblemABC.h"

float CommunicationMinValue(float Value);
int FLDCorrectForImpulses(int field, LevelHierarchyEntry *LevelArray[],
			  int level, Star *AllStars, FLOAT FLDTime, 
			  float dtFLD);

int RadiativeTransferCallFLD(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *AllStars, 
			     ImplicitProblemABC *ImplicitSolver)

{

  LevelHierarchyEntry *Temp;
  int l, CallLevel;
  float dtLocal, dtGrid;
  FLOAT LevelTime;
  Star *cstar;

  if (RadiativeTransferFLD == FALSE)
    return SUCCESS;

  /* Define which field we're calculating.  We should really have some
     logic to determine which one we're calculating. */
  int FieldToCalculate = kdissH2I;

  /* 
     Determine whether we're the coarsest level that is
     >=FLDCallOnLevel (L*) on the same time as level L*.  If so, we
     call the solver.  Below the Xs indicate where we call the solver
     if L* = 2.

     Level 0: X-----------------------X
     Level 1: |-----------X-----------|
     Level 2: |-----X-----|-----X-----|

  */

  for (l = level; l <= RadiativeTransferFLDCallOnLevel; l++) {
    LevelTime = LevelArray[l]->GridData->ReturnTime();
    if (LevelTime >= MetaData->FLDTime) {
      CallLevel = l;
      break;
    }
  }


  // Mode 1: FLD solver used to propagate free-streaming LW radiation 
  //         for ray-tracing solver
  if (RadiativeTransferFLD == 1) {
    if (level == CallLevel) {

      /* If this level is less than the FLD timestepping level, we have
	 to calculate the timestep on level L*.  Otherwise, it's already
	 stored in dtFixed. */

      if (level < RadiativeTransferFLDCallOnLevel) {
	dtLocal = huge_number;
	for (Temp = LevelArray[RadiativeTransferFLDCallOnLevel];
	     Temp; Temp = Temp->NextGridThisLevel) {
	  dtGrid = Temp->GridData->ComputeTimeStep();
	  dtLocal = min(dtLocal, dtGrid);
	}
	MetaData->dtFLD = CommunicationMinValue(dtLocal);
      } else
	MetaData->dtFLD = LevelArray[RadiativeTransferFLDCallOnLevel]->
	  GridData->ReturnTimeStep();
      
      if (debug)
	fprintf(stdout, "CallFLD: level %"ISYM", dtFLD = %g\n",
		level, MetaData->dtFLD);
      
      /* Construct emissivity field on each root grid from the star
	 particles (objects) */
      
      for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->CreateEmissivityLW(AllStars, MetaData->FLDTime, 
					   MetaData->dtFLD);
      
      /* Call FLD solver */
      
      for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
	ImplicitSolver->Evolve(Temp->GridHierarchyEntry, MetaData->dtFLD);
      
      /* Delete emissivity fields */
      
      for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->DeleteEmissivity();
      
      /* For subgrids, we need to interpolate the radiation field from
	 its parent */
      
      HierarchyEntry *Parent;
      for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++)
	for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
	  Parent = Temp->GridHierarchyEntry->ParentGrid;
	  Temp->GridData->InterpolateRadiationFromParent(Parent->GridData, 
							 FieldToCalculate);
	}
      
      MetaData->FLDTime += MetaData->dtFLD;
      
    } // ENDIF level == CallOnLevel
    
    /* On other levels, we have to check for newly-created stars or
       sources between calls to the FLD solver and add a temporary 1/r^2
       radiation field. */
    
    else {
      
      FLDCorrectForImpulses(FieldToCalculate, LevelArray, level, AllStars,
			    MetaData->FLDTime, MetaData->dtFLD);
      
    } // ENDELSE level == CallOnLevel
  } // ENDIF RadiativeTransfer == 1
  

  // For FLD-only solves, we currently call the solver on only the 
  // relevant level (and ignore interpolations to parents).  This may 
  // only work for unigrid runs, or on statically-nested runs where the 
  // FLD solver is called on a single level only.  It should be
  // extended to the general AMR case when the time comes.  However, we 
  // must take care since the HYPRE solver requires that all processes 
  // in the communicator be involved in each solve, and for now we only 
  // have MPI_COMM_WORLD (so not all processes interact on all implicit 
  // "solves" on subgrids).
  else {
    if (level == CallLevel) {
      // use the hydro time step size from this grid level
      if (level < RadiativeTransferFLDCallOnLevel) {
	dtLocal = huge_number;
	for (Temp = LevelArray[RadiativeTransferFLDCallOnLevel];
	     Temp; Temp = Temp->NextGridThisLevel) {
	  dtGrid = Temp->GridData->ComputeTimeStep();
	  dtLocal = min(dtLocal, dtGrid);
	}
	MetaData->dtFLD = CommunicationMinValue(dtLocal);
      } else
	MetaData->dtFLD = LevelArray[RadiativeTransferFLDCallOnLevel]->
	  GridData->ReturnTimeStep();
      
      if (debug)
	fprintf(stdout, "CallFLD: level %"ISYM", dtFLD = %g\n",
		level, MetaData->dtFLD);

      
      // Construct emissivity field on grid from star particles (Geoffrey?)

      
      // Call FLD solver
      for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
	ImplicitSolver->Evolve(Temp->GridHierarchyEntry, MetaData->dtFLD);

      
      // If we made an emissivity field earlier, delete it here

      
      // Update the FLD time
      MetaData->FLDTime += MetaData->dtFLD;
      
    } // ENDIF level == CallOnLevel
  } // ENDELSE RadiativeTransfer == 1
    
  return SUCCESS;

}
