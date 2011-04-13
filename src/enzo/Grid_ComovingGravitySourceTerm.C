/***********************************************************************
/
/  GRID CLASS (ADJUST THE GRAVITY SOURCE TERM FOR COMOVING COORDINATES)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int grid::ComovingGravitySourceTerm()
{
 
  /* Return is this is not the right processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* This assumes that all the equations have units as defined by
     CosmologyUnits.  This means that GravitationalConstant must be 1. */
 
  if (GravitationalConstant != 1) {
    ENZO_FAIL("GravitationalConstant must be 1!.\n");
  }
 
  /* Set AverageDensity (held in global_data.h). */
 
  int i, j, k, gmf_index;
  float AverageDensity = 1.;

  if (ProblemType == 50 || ProblemType == 60 || ProblemType == 61) //AK
    AverageDensity = 0.0;
 
  /* Loop over the field, subracting off the mean field. */
 
  for (k = 0; k < GravitatingMassFieldDimension[2]; k++)
    for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
      gmf_index = (k*GravitatingMassFieldDimension[1] + j)*
	GravitatingMassFieldDimension[0];
      for (i = 0; i < GravitatingMassFieldDimension[0]; i++, gmf_index++) {
	GravitatingMassField[gmf_index] =
	  GravitatingMassField[gmf_index] - AverageDensity;
      }
    }
 
  return SUCCESS;
}
