/***********************************************************************
/
/  PROTOSUBGRID CLASS (RETURNS THE NTH LONGEST GRID DIMENSION)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
int ProtoSubgrid::ReturnNthLongestDimension(int n)
{
 
  int dim, TempInt[MAX_DIMENSION], Order[MAX_DIMENSION], i, swap;
 
  /* Return -1 if N is too large (N is zero-based). */
 
  if (n >= GridRank)
    return -1;
 
  /* Copy dimensions into temp and sort. */
 
  for (dim = 0; dim < GridRank; dim++) {
    TempInt[dim] = GridDimension[dim];
    Order[dim] = dim;
  }
 
  /* Sort. */
 
  for (i = 0; i < GridRank-1; i++)
    for (dim = 0; dim < GridRank-1; dim++)
      if (TempInt[dim] < TempInt[dim+1]) {
	swap = TempInt[dim+1];
	TempInt[dim+1] = TempInt[dim];
	TempInt[dim] = swap;
	swap = Order[dim+1];
	Order[dim+1] = Order[dim];
	Order[dim] = swap;
      }
 
  /* Return the Nth longest dim. */
 
  return Order[n];
}
