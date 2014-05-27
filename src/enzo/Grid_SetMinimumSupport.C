/***********************************************************************
/
/  GRID CLASS (SET THE ENERGY TO PROVIDE MINIMAL PRESSURE SUPPORT)
/
/  written by: Greg Bryan
/  date:       November, 1998
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
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::SetMinimumSupport(float &MinimumSupportEnergyCoefficient)
{

  if (NumberOfBaryonFields > 0) {
 
    const float pi = 3.14159;
 
    /* Compute cosmology factors. */
 
    FLOAT a = 1, dadt;
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
	ENZO_FAIL("Error in CosmologyComputeExpansionFactor.\n");
      }
    float CosmoFactor = 1.0/a;
 
    /* Calculate scaling factor for minimum pressure */
 
    MinimumSupportEnergyCoefficient =
      GravitationalConstant/(4.0*pi) / (pi * (Gamma*(Gamma-1.0))) *
                  CosmoFactor * MinimumPressureSupportParameter *
                  CellWidth[0][0] * CellWidth[0][0];
 
  } // end: if (NumberOfBaryonFields > 0)

 
  return SUCCESS;
}
