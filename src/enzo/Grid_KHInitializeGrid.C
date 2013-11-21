/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR KH PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, December 2004.  
/  modified2:  Gregg Dobrowalski, Feb 2005.  
/  modified3:  Alexei Kritsuk, April 2005. v_y perturbations + more parameters.
/  modified4:  Cameron Hummels, November 2013. Convergent ICs
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

int grid::KHInitializeGrid(float KHInnerDensity, 
			   float KHOuterDensity,
			   float KHInnerInternalEnergy,
			   float KHOuterInternalEnergy,
			   float KHPerturbationAmplitude,
			   float KHInnerVx, float KHOuterVx,
			   float KHInnerPressure,
			   float KHOuterPressure,
               float KHConvergentICs,
               float KHRampWidth)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  if (KHConvergentICs == 0.0) {   // use old method with divergent behavior
                                  // which sets up a step function between
                                  // the two fluids and perturbs with random
                                  // fluctutations

    /* set fields in "inner" region: .25 < y < .75 */

    int index, jndex, i;
    for (i = 0; i < size; i++) {
      index = i % GridDimension[0];
      jndex = (i-index)/GridDimension[0];

      BaryonField[3][i] = KHPerturbationAmplitude * 
                          ((float)rand()/(float)(RAND_MAX) - 0.5); //AK
      BaryonField[1][i] = POW(BaryonField[3][i],2) / 2.0;
      
      if (jndex >= GridDimension[1]/4 && jndex < 3*GridDimension[1]/4) { //AK
        BaryonField[0][i] = KHInnerDensity;
        BaryonField[2][i] = KHInnerVx + 
	                        KHPerturbationAmplitude * 
	                        ((float)rand()/(float)(RAND_MAX) - 0.5);
        BaryonField[1][i] += KHInnerInternalEnergy + 
	                         POW(BaryonField[2][i],2) / 2.0;
      }
      else {
        BaryonField[0][i] = KHOuterDensity;
        BaryonField[2][i] = KHOuterVx + 
	                        KHPerturbationAmplitude * 
                             ((float)rand()/(float)(RAND_MAX) - 0.5);
        BaryonField[1][i] += KHOuterInternalEnergy + 
	                         POW(BaryonField[2][i],2) / 2.0;
      }
    }
  }
  else {    // use convergent initial conditions, which create a continuously
            // changing ``ramp`` region between the two fluids in density/vel
            // space.  then perturb the fluids with a sinusoidal perturbation.
            // all resulting simulations should converge as you increase 
            // resolution
    int index, jndex, n;
    FLOAT x, y, KHRampWidth2;
    KHRampWidth2 = KHRampWidth*KHRampWidth;
    for (int j = 0; j < GridDimension[1]; j++)
      for (int i = 0; i < GridDimension[0]; i++, n++) {

     /* Compute physical position of each cell */

     x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
     y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];

//      fprintf(stderr, "CellLeftEdge = (%"FSYM" %"FSYM")\n", CellLeftEdge[0][i], CellLeftEdge[1][j]);
//    fprintf(stderr, "x = %"FSYM" y = %"FSYM"\n", x, y);

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
      } else {
      // set density and velocity fields in "outer region" 
        BaryonField[0][n] = KHOuterDensity;
        BaryonField[2][n] = KHOuterVx;
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
  }
      
  if (debug) {
    fprintf(stderr, "GKHIG: BF[2][0]:%"FSYM" BF[2][size/2]: %"FSYM"\n", 
	    BaryonField[2][0],
	    BaryonField[2][size/2]);
  
    fprintf(stderr, "GKHIG: BF[3][0]:%"FSYM" BF[3][size/2]: %"FSYM"\n", 
	    BaryonField[3][0],
	    BaryonField[3][size/2]);
  }
  return SUCCESS;
}
