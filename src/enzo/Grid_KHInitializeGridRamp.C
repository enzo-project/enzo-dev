/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR KH PROBLEM)
/
/  written by: Cameron Hummels
/  date:       November, 2013
/
/  PURPOSE: Sets the field variables in the domain. Adds ramp between 
/           the two fluids (in density and x-velocity) so transition is 
/           continuous.  Also perturb y-velocity field sinusoidally to 
/           seed KHI. These simulations converge in behavior with increasing
/           resolution unlike the non-ramp method (i.e. KHInitializeGrid).
/           These ICs are based on the paper by:
/           Robertson, Kravtsov, Gnedin, Abel & Rudd 2010
/           N.B. The paper has a typo in the ICs ramp equation, got original
/           from Brant Robertson for this prescription.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::KHInitializeGridRamp(float KHInnerDensity, 
                               float KHOuterDensity,
                               float KHInnerInternalEnergy,
                               float KHOuterInternalEnergy,
                               float KHPerturbationAmplitude,
                               float KHInnerVx, float KHOuterVx,
                               float KHInnerPressure,
                               float KHOuterPressure,
                               float KHRampWidth)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this-IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                      Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  int index, jndex, i, j, k, n = 0;
  FLOAT x, y, KHRampWidth2, DensityDiff, VelocityDiff;
  KHRampWidth2 = KHRampWidth*KHRampWidth;
  DensityDiff = (KHInnerDensity - KHOuterDensity);
  VelocityDiff = (KHInnerVx - KHOuterVx);
  for (k = 0; k < GridDimension[2]; k++)
  for (j = 0; j < GridDimension[1]; j++)
  for (i = 0; i< GridDimension[0]; i++, n++) {

    /* Define coordinates based on absolute frame */
    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];

    /* everywhere set the sinusoidal perturbation in v_y */
    BaryonField[Vel2Num][n] = KHPerturbationAmplitude * sin(4.0 * M_PI * x);

    /* set density and velocity fields in "inner region" */
    if (y > 0.25 && y < 0.75) { 
      BaryonField[DensNum][n]  = KHInnerDensity;
      BaryonField[Vel1Num][n]  = KHInnerVx;

    /* modify the top half of the inner fluid to account for ramp */
      if (y >= 0.5) {
        BaryonField[DensNum][n] -= DensityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[Vel1Num][n] -= VelocityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {

    /* modify the bottom half of the inner fluid to account for ramp */
        BaryonField[DensNum][n] -= DensityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[Vel1Num][n] -= VelocityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }

    /* set the energy in the "inner region" to be thermal + kinetic */
      BaryonField[TENum][n] = KHInnerPressure/((Gamma - 1.0)*BaryonField[DensNum][n]) + 
                       (POW(BaryonField[Vel1Num][n],2) +
                        POW(BaryonField[Vel2Num][n],2)) / 2.0;

    /* set density and velocity fields in "outer region" */
    } else {
      BaryonField[DensNum][n]  = KHOuterDensity;
      BaryonField[Vel1Num][n]  = KHOuterVx;

    /* modify the top half of the outer fluid to account for ramp */
      if (y >= 0.5) {
        BaryonField[DensNum][n] += DensityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[Vel1Num][n] += VelocityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {

    /* modify the bottom half of the outer fluid to account for ramp */
        BaryonField[DensNum][n] += DensityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[Vel1Num][n] += VelocityDiff * exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }

    /* set the energy in the "outer region" to be thermal + kinetic */
      BaryonField[TENum][n] = KHOuterPressure/((Gamma - 1.0)*BaryonField[DensNum][n]) + 
                       (POW(BaryonField[Vel1Num][n],2) +
                        POW(BaryonField[Vel2Num][n],2)) / 2.0;
    }
  }
  return SUCCESS;
}
