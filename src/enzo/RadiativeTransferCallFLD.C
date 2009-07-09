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

int RadiativeTransferCallFLD(LevelHierarchyEntry *LevelArray[], int level,
			     Star *AllStars, ImplicitProblemABC *ImplicitSolver)

{

  LevelHierarchyEntry *Temp;
  FLOAT FLDTime;
  float dt;

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

  int l, CallLevel;
  FLOAT LevelTime;
  FLDTime = LevelArray[RadiativeTransferFLDCallOnLevel]->GridData->
    ReturnTime();
  for (l = 0; l <= RadiativeTransferFLDCallOnLevel; l++) {
    LevelTime = LevelArray[l]->GridData->ReturnTime();
    if (LevelTime == FLDTime) {
      CallLevel = l;
      break;
    }
  }

  if (level == CallLevel) {

    dt = LevelArray[level]->GridData->ReturnTimeStep();

    /* Construct emissivity field on each root grid from the star
       particles (objects) */

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->CreateEmissivityLW(AllStars, FLDTime, dt);

    /* Call FLD solver */

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      ImplicitSolver->Evolve(Temp->GridHierarchyEntry, dt);

    /* Delete emissivity fields */

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->DeleteEmissivity();

    /* For subgrids, we need to interpolate the radiation field from
       its parent */

    HierarchyEntry *Parent;
    for (level = 1; level < MAX_DEPTH_OF_HIERARCHY; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
	Parent = Temp->GridHierarchyEntry->ParentGrid;
	Temp->GridData->InterpolateRadiation(Parent->GridData, FieldToCalculate);
      }

  } // ENDIF level == CallOnLevel

  return SUCCESS;

}
