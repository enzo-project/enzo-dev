/***********************************************************************
/
/  GRID CLASS (WRAP THE GRACKLE CHEMISTRY SOLVER)
/
/  written by: Britton Smith
/  date:       April, 2013
/  modified1:
/
/  PURPOSE: Solve chemistry and cooling with grackle.
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
 
int grid::GrackleWrapper()
{

  if (!grackle_chemistry.use_chemistry)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  LCAPERF_START("grid_GrackleWrapper");

  LCAPERF_STOP("grid_GrackleWrapper");
  return SUCCESS;
}
