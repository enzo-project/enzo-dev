/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR KH PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, December 2004.  
/  modified2:  Gregg Dobrowalski, Feb 2005.  
/  modified3:  Alexei Kritsuk, April 2005. v_y perturbations + more parameters.
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
			   float KHInnerInternalEnergy,
			   float KHOuterInternalEnergy,
			   float KHPerturbationAmplitude,
			   float KHInnerVx, float KHOuterVx)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in "inner" region: .25 < y < .75 */

  int index, jndex, i;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    BaryonField[3][i] = KHPerturbationAmplitude * 
                        ((float)rand()/(float)(RAND_MAX) - 0.5); //AK
    BaryonField[1][i] = POW(BaryonField[3][i],2) / 2.0;
    
    if (jndex >= GridDimension[1]/4 && jndex < 3*GridDimension[1]/4) { //AK
      BaryonField[0][i]  = KHInnerDensity;
      BaryonField[2][i]  = KHInnerVx + 
	                   KHPerturbationAmplitude * 
	                   ((float)rand()/(float)(RAND_MAX) - 0.5);
      BaryonField[1][i] += KHInnerInternalEnergy + 
	                   POW(BaryonField[2][i],2) / 2.0;
    }
    else {
      BaryonField[2][i]  = KHOuterVx + 
	                   KHPerturbationAmplitude * 
                           ((float)rand()/(float)(RAND_MAX) - 0.5);
      BaryonField[1][i] += KHOuterInternalEnergy + 
	                   POW(BaryonField[2][i],2) / 2.0;
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

















































