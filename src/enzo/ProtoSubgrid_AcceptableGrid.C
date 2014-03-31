/***********************************************************************
/
/  PROTOSUBGRID CLASS (CHECK IF A SUBGRID NEEDS TO BE SUBDIVIDED)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:  Robert Harkness
/  date:       November, 2005
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def" 
 
/* function prototypes */

/*
   These were definitions in the original version but have been
   superseded by parameters for tuning in large-scale AMR runs

   MINIMUM_SIDE_LENGTH 4
   MAXIMUM_SIZE 2000
*/
 
int ProtoSubgrid::AcceptableSubgrid()
{
 
  /* If NumberFlagged hasn't been computed yet, then compute it. */
 
  if (NumberFlagged == INT_UNDEFINED) {
    this->ComputeSignature(0);
    NumberFlagged = 0;
    for (int i = 0; i < GridDimension[0]; i++)
      NumberFlagged += Signature[0][i];
  }

  /* Compute size and efficiency. */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
    if( GridDimension[dim] >= 0.5*MAX_ANY_SINGLE_DIRECTION)
        return FALSE;
  }
 
  float efficiency = float(NumberFlagged)/float(size);

  /* For size restrictions, use the number of child cells */
 
  //size *= 8;

  /* Subgrid is acceptable if it is efficient enough or small enough,
     but if we are using multiple processors, then make sure it is
     not too big (for good load balancing). */
 
  int DimLong = this->ReturnNthLongestDimension(0);
  int DimShort = this->ReturnNthLongestDimension(GridRank-1);

  float ratio = float(GridDimension[DimLong])/float(GridDimension[DimShort]);

  //if (ratio > 3.0)
  //  return FALSE;

  if (size <= POW(float(MinimumSubgridEdge), GridRank)){
    if (ratio > 3.0)
      printf("Size acceptable; ratio is %g\n",ratio);
    return TRUE;
  }
 
  if (size > MaximumSubgridSize && NumberOfProcessors > 1)
    return FALSE;
 
  if (efficiency > MinimumEfficiency){
    if (ratio > 3.0)
      printf("Efficiency acceptable; ratio is %g; efficiency is %g > %g\n",ratio,efficiency,MinimumEfficiency);
    return TRUE;
  }
 
  /* Not acceptable yet -- keep refining. */
 
  return FALSE;
}
