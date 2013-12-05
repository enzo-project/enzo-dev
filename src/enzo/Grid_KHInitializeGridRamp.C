/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR KH PROBLEM)
/
/  written by: Cameron Hummels
/  date:       November, 2013
/
/  PURPOSE: Sets the field variables in the domain.
/          
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
                               float KHRampWidth,
                               int   level)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int index, jndex, i, j, k, n = 0;
  FLOAT x, y, KHRampWidth2;
  KHRampWidth2 = KHRampWidth*KHRampWidth;
  for (k = 0; k < GridDimension[2]; k++)
  for (j = 0; j < GridDimension[1]; j++)
  for (i = 0; i< GridDimension[0]; i++, n++) {
/*  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    x = (FLOAT)(index - GridStartIndex[0]) / 
            (FLOAT)(GridEndIndex[0] - GridStartIndex[0]+1);
    y = (FLOAT)(jndex - GridStartIndex[1]) / 
            (FLOAT)(GridEndIndex[1] - GridStartIndex[1]+1);
*/
    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
//    fprintf(stderr, "i = %"ISYM" j = %"ISYM", n = %"ISYM" x = %"FSYM" y = %"FSYM"\n", i, j, n, x, y);

    // everywhere set the sinusoidal perturbation in v_y
    BaryonField[3][n] = KHPerturbationAmplitude * sin(4.0 * M_PI * x);

    // set density and velocity fields in "inner region"
    if (y > 0.25 && y < 0.75) { 
      BaryonField[0][n]  = KHInnerDensity;
      BaryonField[2][n]  = KHInnerVx;
    // modify the top half of the inner fluid (quadrant 3) to account for ramp
      if (y >= 0.5) {
        BaryonField[0][n] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][n] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {
    // modify the bottom half of the inner fluid (quadrant 2) to account for ramp
        BaryonField[0][n] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][n] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }
    // set the energy in the "inner region" 
      BaryonField[1][n] = KHInnerPressure/((Gamma - 1.0)*BaryonField[0][n]) + 
                       (POW(BaryonField[2][n],2) +
                        POW(BaryonField[3][n],2)) / 2.0;

    // set density and velocity fields in "outer region"
    } else {
      BaryonField[0][n]  = KHOuterDensity;
      BaryonField[2][n]  = KHOuterVx;
    // modify the top half of the outer fluid (quadrant 4) to account for ramp
      if (y >= 0.5) {
        BaryonField[0][n] += exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][n] += exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {
    // modify the bottom half of the outer fluid (quadrant 1) to account for ramp
        BaryonField[0][n] += exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][n] += exp( (-0.5/KHRampWidth2)*
                                  pow(y-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }
    // set the energy in the "outer region"
      BaryonField[1][n] = KHOuterPressure/((Gamma - 1.0)*BaryonField[0][n]) + 
                       (POW(BaryonField[2][n],2) +
                        POW(BaryonField[3][n],2)) / 2.0;
    }
  }
  return SUCCESS;
}
