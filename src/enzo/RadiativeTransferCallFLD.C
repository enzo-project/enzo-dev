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

  /* Return if not root level or not needed */

  if (RadiativeTransferFLD == FALSE)
    return SUCCESS;

  LevelHierarchyEntry *Temp;

  /* Construct emissivity field on each root grid from the star
     particles (objects) */
  
  if (level == 0) {

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->CreateEmissivityLW(AllStars);

    /* Call FLD solver */

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      ImplicitSolver->Evolve(Temp->GridHierarchyEntry);

    /* Delete emissivity fields */

    for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->DeleteEmissivity();

  } // ENDIF level == 0

  /* For subgrids, we need to interpolate the radiation field from its
     parent */

  else {

    

  } // ENDELSE level == 0
  
  return SUCCESS;

}
