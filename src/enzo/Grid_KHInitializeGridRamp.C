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
               float KHRampWidth)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int index, jndex, i;
  FLOAT x_val, y_val, KHRampWidth2;
  KHRampWidth2 = KHRampWidth*KHRampWidth;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    x_val = (FLOAT)(index - GridStartIndex[0]) / 
            (FLOAT)(GridEndIndex[0] - GridStartIndex[0]+1);
    y_val = (FLOAT)(jndex - GridStartIndex[1]) / 
            (FLOAT)(GridEndIndex[1] - GridStartIndex[1]+1);

    // everywhere set the sinusoidal perturbation in v_y
    BaryonField[3][i] = KHPerturbationAmplitude * sin(4.0 * M_PI * x_val);

    // set density and velocity fields in "inner region"
    if (y_val > 0.25 && y_val < 0.75) { 
      BaryonField[0][i]  = KHInnerDensity;
      BaryonField[2][i]  = KHInnerVx;
    // modify the top half of the inner fluid (quadrant 3) to account for ramp
      if (y_val >= 0.5) {
        BaryonField[0][i] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][i] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.75 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {
    // modify the bottom half of the inner fluid (quadrant 2) to account for ramp
        BaryonField[0][i] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][i] -= exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.25 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }
    // set the energy in the "inner region" 
      BaryonField[1][i] = KHInnerPressure/((Gamma - 1.0)*BaryonField[0][i]) + 
                       (POW(BaryonField[2][i],2) +
                        POW(BaryonField[3][i],2)) / 2.0;
    } else {
    // density and velocity fields in "outer region" already set in KHInit...C
    // modify the top half of the outer fluid (quadrant 4) to account for ramp
      if (y_val >= 0.5) {
        BaryonField[0][i] += exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][i] += exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.75 + 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      } else {
    // modify the bottom half of the outer fluid (quadrant 1) to account for ramp
        BaryonField[0][i] += exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
        BaryonField[2][i] += exp( (-0.5/KHRampWidth2)*
                                  pow(y_val-0.25 - 
                                      sqrt(-2.0*KHRampWidth2*log(0.5)),2));
      }
    // set the energy in the "outer region"
      BaryonField[1][i] = KHOuterPressure/((Gamma - 1.0)*BaryonField[0][i]) + 
                       (POW(BaryonField[2][i],2) +
                        POW(BaryonField[3][i],2)) / 2.0;
    }
  }
  return SUCCESS;
}
