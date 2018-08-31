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

#include "preincludes.h"
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

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {
    grackle_data->radiative_transfer_intermediate_step = FALSE;
    if (this->GrackleWrapper() == FAIL) {
      ENZO_FAIL("Error in GrackleWrapper.\n");
    }
    return SUCCESS;
  }
#endif

  if (MultiSpecies && RadiativeCooling ) {
    int RTCoupledSolverIntermediateStep = FALSE;
    this->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
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
