/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR KH PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, December 2004.  
/  modified2:  Gregg Dobrowalski, Feb 2005.  
/  modified3:  Alexei Kritsuk, April 2005. v_y perturbations + more parameters.
/  modified4:  Cameron Hummels, December 2013. mersennes twister RNG
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

/* including mersennes twister random number generator; allows
   reproducible behavior with same seed */

void mt_init(unsigned_int seed);

unsigned_long_int mt_random();

int grid::KHInitializeGrid(float KHInnerDensity, 
                           float KHInnerInternalEnergy,
                           float KHOuterInternalEnergy,
                           float KHPerturbationAmplitude,
                           float KHInnerVx, float KHOuterVx,
                           float KHInnerPressure,
                           float KHOuterPressure,
                           int   KHRandomSeed)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  float rand_x, rand_y;
  unsigned_long_int random_int;

  /* initialize random number generator from seed */
  mt_init(((unsigned_int) KHRandomSeed)); 

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

    /* set fields in "inner" region: .25 < y < .75 */

  int index, jndex, i;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];

    rand_x = (float)(mt_random()%32768)/(32768.0) - 0.5;
    rand_y = (float)(mt_random()%32768)/(32768.0) - 0.5;

    BaryonField[3][i] = KHPerturbationAmplitude * rand_y;
    BaryonField[1][i] = POW(BaryonField[3][i],2) / 2.0;
      
    if (jndex >= GridDimension[1]/4 && jndex < 3*GridDimension[1]/4) { //AK
      BaryonField[0][i]  = KHInnerDensity;
      BaryonField[2][i]  = KHInnerVx + KHPerturbationAmplitude * rand_x;
      BaryonField[1][i] += KHInnerInternalEnergy + 
                           POW(BaryonField[2][i],2) / 2.0;
    }
    else {
      BaryonField[2][i]  = KHOuterVx + KHPerturbationAmplitude * rand_x;
      BaryonField[1][i] += KHOuterInternalEnergy + 
                           POW(BaryonField[2][i],2) / 2.0;
    }
  }
      
  if (debug) {
    fprintf(stderr, "GKHIG: BF[2][0]:%"FSYM" BF[2][size/2]: %"FSYM"\n", 
            BaryonField[2][0], BaryonField[2][size/2]);

    fprintf(stderr, "GKHIG: BF[3][0]:%"FSYM" BF[3][size/2]: %"FSYM"\n", 
            BaryonField[3][0], BaryonField[3][size/2]);
  }
  return SUCCESS;
}
