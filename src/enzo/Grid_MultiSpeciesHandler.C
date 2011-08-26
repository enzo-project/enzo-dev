/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING COOLING/CHEMISTRY)
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: Move logic for chemistry/cooling module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

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
 
int grid::MultiSpeciesHandler()
{
  if ((!MultiSpecies) && (!RadiativeCooling)) return SUCCESS; 
  if (GadgetEquilibriumCooling != 0) return SUCCESS;

  LCAPERF_START("grid_MultiSpeciesHandler");

  if (MultiSpecies && RadiativeCooling ) {
    if((MultiSpecies == 3) && (PrimordialChemistrySolver > 0)) {
      if (PrimordialChemistrySolver == 1) {
	this->SolveHighDensityPrimordialChemistry();
      } else if (PrimordialChemistrySolver == 2) {
#ifdef USE_CVODE
	this->SolvePrimordialChemistryCVODE();
#else
	ENZO_FAIL("You have asked for the CVODE solver but CVODE not enabled!");
#endif
      }
    } // ENDIF PrimordialChemistrySolver
    else {
      int RTCoupledSolverIntermediateStep = FALSE;
      this->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
    }
  } else {
    if (MultiSpecies)
      this->SolveRateEquations();
    if (RadiativeCooling)
      this->SolveRadiativeCooling();
  }

  if (ProblemType == 62)
    this->CoolingTestResetEnergies();

  LCAPERF_STOP("grid_MultiSpeciesHandler");
  return SUCCESS;
}
