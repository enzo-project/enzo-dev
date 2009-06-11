/***********************************************************************
/
/  PROTOSUBGRID CLASS (COMPUTE SECOND DERIVATIVE OF SIGNATURE AND FIND
/                      THE STRONGEST ZERO CROSSING)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int ProtoSubgrid::ComputeSecondDerivative(int CheckDim,
					  int &ZeroCrossStrength,
					  int GridEnds[2][2])
{
  /* Initialize. */
 
  ZeroCrossStrength = 0;
  int i, Strength, Center = (GridDimension[CheckDim]-1)/2;
  int ZeroCrossPosition = Center;
  GridEnds[0][0] = StartIndex[CheckDim];
  GridEnds[1][1] = EndIndex[CheckDim];
 
  /* If the dimension is < 4, then we can't compute the 2nd derivative, so
     return the zero strength crossing with a break in the middle. */
 
  if (GridDimension[CheckDim] < 4) {
    GridEnds[0][1] = StartIndex[CheckDim] + Center;
    GridEnds[1][0] = min(GridEnds[0][1]+1, EndIndex[CheckDim]);
    return SUCCESS;
  }
 
  /* Compute the signature along CheckDim, if necessary. */
 
  this->ComputeSignature(CheckDim);
 
  /* Allocate a new arrays to hold the derivative values. */
 
  int *SecondDerivative = new int[GridDimension[CheckDim]];
 
  /* Compute d2q/dx2. */
 
  for (i = 1; i < GridDimension[CheckDim]-1; i++)
    SecondDerivative[i] =   Signature[CheckDim][i-1] -
                          2*Signature[CheckDim][i  ] +
                            Signature[CheckDim][i+1] ;
 
  /* Find strongest zero crossing with ties won by the position closest to the
     center. */
 
  for (i = 1; i < GridDimension[CheckDim]-2; i++)
    if (SecondDerivative[i]*SecondDerivative[i+1] <= 0) {
      Strength = ABS(SecondDerivative[i] - SecondDerivative[i+1]);
      if (Strength > ZeroCrossStrength ||
	  (Strength == ZeroCrossStrength && ABS(Center - i) <
	                                    ABS(Center - ZeroCrossPosition))) {
	ZeroCrossStrength = Strength;
	ZeroCrossPosition = i;
      }
    }
 
  /* Now, set the rest of the grid start and end indicies. */
 
  GridEnds[0][1] = GridEnds[0][0] + ZeroCrossPosition;
  GridEnds[1][0] = GridEnds[0][0] + ZeroCrossPosition + 1;
 
  /* Clean up */
 
  delete [] SecondDerivative;
 
  return SUCCESS;
}
