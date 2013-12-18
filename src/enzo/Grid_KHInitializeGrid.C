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

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this-IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                      Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* set fields in "inner" region: .25 < y < .75 */

  int index, jndex, i;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];

    rand_x = (float)(mt_random()%32768)/(32768.0) - 0.5;
    rand_y = (float)(mt_random()%32768)/(32768.0) - 0.5;

    BaryonField[Vel2Num][i] = KHPerturbationAmplitude * rand_y;
    BaryonField[TENum][i] = POW(BaryonField[Vel2Num][i],2) / 2.0;
      
    if (jndex >= GridDimension[1]/4 && jndex < 3*GridDimension[1]/4) { //AK
      BaryonField[DensNum][i]  = KHInnerDensity;
      BaryonField[Vel1Num][i]  = KHInnerVx + KHPerturbationAmplitude * rand_x;
      BaryonField[TENum][i] += KHInnerInternalEnergy + 
                               POW(BaryonField[Vel1Num][i],2) / 2.0;
    }
    else {
      BaryonField[Vel1Num][i]  = KHOuterVx + KHPerturbationAmplitude * rand_x;
      BaryonField[TENum][i] += KHOuterInternalEnergy + 
                               POW(BaryonField[Vel1Num][i],2) / 2.0;
    }
  }
      
  if (debug) {
    fprintf(stderr, "GKHIG: BF[2][0]:%"FSYM" BF[2][size/2]: %"FSYM"\n", 
            BaryonField[Vel1Num][0], BaryonField[Vel1Num][size/2]);

    fprintf(stderr, "GKHIG: BF[3][0]:%"FSYM" BF[3][size/2]: %"FSYM"\n", 
            BaryonField[Vel2Num][0], BaryonField[Vel2Num][size/2]);
  }
  return SUCCESS;
}
