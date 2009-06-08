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
	fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
	return FAIL;
      }
    float CosmoFactor = 1.0/a;
 
    /* Determine the size of the grids. */
 
    int dim, size = 1, i;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Find the density, gas energy, velocities & total energy
       (where appropriate). */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }
 
    /* Set minimum GE. */
 
    MinimumSupportEnergyCoefficient =
      GravitationalConstant/(4.0*pi) / (pi * (Gamma*(Gamma-1.0))) *
                  CosmoFactor * MinimumPressureSupportParameter *
                  CellWidth[0][0] * CellWidth[0][0];
 
    if (HydroMethod == Zeus_Hydro)
      ;
 
      /* Zeus: TE is really GE (ugh) -- now applied to pressure inside
	   zeus routine (much better and should be applied to ppm as well)
 
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] = max(BaryonField[TENum][i], Coeff*
				    BaryonField[DensNum][i]);   */
    else {
 
      /* PPM: set GE. */
 
      if (DualEnergyFormalism == TRUE) {
	for (i = 0; i < size; i++)
	  BaryonField[GENum][i] = max(BaryonField[GENum][i],
				      MinimumSupportEnergyCoefficient *
				      BaryonField[DensNum][i]);
	if (GridRank != 3) return FAIL;
	for (i = 0; i < size; i++)
	  BaryonField[TENum][i] = max(BaryonField[GENum][i] + 0.5*
		       (BaryonField[Vel1Num][i]*BaryonField[Vel1Num][i] +
		        BaryonField[Vel2Num][i]*BaryonField[Vel2Num][i] +
		        BaryonField[Vel3Num][i]*BaryonField[Vel3Num][i]),
				      BaryonField[TENum][i]);
								
      }
      else {
	fprintf(stderr, "not implemented.\n");
	return FAIL;
      }
    }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  return SUCCESS;
}
