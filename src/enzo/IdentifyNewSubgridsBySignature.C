/***********************************************************************
/
/  RECURSIVELY IDENTIFY NEW SUBGRIDS BY SIGNATURES
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
static int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2];
 
int IdentifyNewSubgridsBySignature(ProtoSubgrid *SubgridList[],
				   int &NumberOfSubgrids)
{
 
  int dim, i, j, NumberOfNewGrids;
  ProtoSubgrid *NewSubgrid, *Subgrid;
 
  /* Loop over all the grids in the queue SubgridList. */

  if ( NumberOfSubgrids > MAX_NUMBER_OF_SUBGRIDS ) {
    fprintf(stderr, "PE %"ISYM" NumberOfSubgrids > MAX_NUMBER_OF_SUBGRIDS in IdentifyNewSubgridsBySignature\n", MyProcessorNumber);
    return FAIL;
  }
 
  int index = 0;

  while (index < NumberOfSubgrids) {
 
    Subgrid = SubgridList[index];
 
    /* Shrink this subgrid (if necessary) to produce the smallest box. */
 
    if (Subgrid->ShrinkToMinimumSize() == FAIL) {
      fprintf(stderr, "Error in ProtoSubgrid->ShrinkToMinimumSize.\n");
      return FAIL;
    }
 
    /* Iterate on this grid until it is acceptable. */
 
    while (Subgrid->AcceptableSubgrid() == FALSE) {
 
      /* Loop over the dimensions (longest to shortest), compute the
	 1D signatures and then look for zeros in them.  */
 
      for (i = 0; i < MAX_DIMENSION; i++) {
 
	/* Find all the zeros in the signature in the ith longest dimension. */
 
	NumberOfNewGrids = 0;

	if ((dim = Subgrid->ReturnNthLongestDimension(i)) < 0)
	  break;

	Subgrid->ComputeSignature(dim);

	if (Subgrid->FindGridsByZeroSignature(dim, NumberOfNewGrids, GridEnds) == FAIL) {
	 fprintf(stderr, "Error in ProtoSubgrid->FindGridsByZeroSignature.\n");
	 return FAIL;
	}
 
	/* Error check. */
 
	if (NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "Increase MAX_NUMBER_OF_SUBGRIDS in IdentifyNewSubgridsBySignature.\n");
	  return FAIL;
	}
 
	/* If there are any new grids created this way, then make them and
	   break out of the loop (note: 1 new grid means no change). */
 
	if (NumberOfNewGrids > 1) {
 
	  /*	  if (debug)
	    printf("Breaking by simple zero. new grids[%"ISYM"]=%"ISYM" break=%"ISYM"\n",
		   dim, NumberOfNewGrids, GridEnds[j][1]); */
 
	  for (j = 0; j < NumberOfNewGrids; j++) {
	    NewSubgrid = new ProtoSubgrid;
	    Subgrid->CopyToNewSubgrid(dim, GridEnds[j][0], GridEnds[j][1],
				      NewSubgrid);
	    if (j == 0)
	      SubgridList[index] = NewSubgrid;
	    else
	      SubgridList[NumberOfSubgrids++] = NewSubgrid;
	  }
 
	  break; // break out of the loop over dimensions
	}
	
      } // end: for (i = 0; i < MAX_DIMENSION; i++)
 
      /* If we couldn't partition by simple zeros then try the 2nd deriv. */
 
      if (NumberOfNewGrids <= 1) {
 
	/* Now Compute the zero crossings in the second derivaties of all the
	   signatures. */
 
	int MaxZeroCrossingStrength = -1, StrongestDim = -1, TempInt;
	
	//	for (i = 0; i < MAX_DIMENSION; i++) {
	for (i = 0; i < 1; i++) {
	
	  if ((dim = Subgrid->ReturnNthLongestDimension(i)) < 0)
	    break;
	
	  if (Subgrid->ComputeSecondDerivative(dim, TempInt,
					       &GridEnds[dim*2]) == FAIL) {
	   fprintf(stderr,"Error in ProtoSubgrid->ComputeSecondDerivative.\n");
	   return FAIL;
	  }
	
	  if (TempInt > MaxZeroCrossingStrength) {
	    StrongestDim = dim;
	    MaxZeroCrossingStrength = TempInt;
	  }
	
	} // end: for (i = 0; i < MAX_DIMENSION; i++)
	
	/* Error check. */
	
	if (StrongestDim < 0) {
	  fprintf(stderr, "Error in IdentifyNewSubgridsBySignature.\n");
	  return FAIL;
	}
	
	/* Create new subgrids (two). */
 
	SubgridList[index] = new ProtoSubgrid;
	SubgridList[NumberOfSubgrids++] = new ProtoSubgrid;
	Subgrid->CopyToNewSubgrid(StrongestDim, GridEnds[StrongestDim*2][0],
				  GridEnds[StrongestDim*2][1],
				  SubgridList[index]);
	Subgrid->CopyToNewSubgrid(StrongestDim, GridEnds[StrongestDim*2+1][0],
				  GridEnds[StrongestDim*2+1][1],
				  SubgridList[NumberOfSubgrids-1]);

	/*	if (debug)
	  printf("Breaking by zero-crossing. dim=%"ISYM"  break=%"ISYM"-%"ISYM"/%"ISYM"-%"ISYM"\n",
		 StrongestDim,
		 GridEnds[StrongestDim*2][0], GridEnds[StrongestDim*2][1],
		 GridEnds[StrongestDim*2+1][0], GridEnds[StrongestDim*2+1][1]);
		 */

      }
 
      /* Delete the old subgrid and set Subgrid to the (first) new grid. */
 
      delete Subgrid;
      Subgrid = SubgridList[index];
 
      /* Shrink this subgrid (if necessary) to produce the smallest box. */
 
      if (Subgrid->ShrinkToMinimumSize() == FAIL) {
	fprintf(stderr, "Error in ProtoSubgrid->ShrinkToMinimumSize.\n");
	return FAIL;
      }
 
    } // end: while (Subgrid->AcceptableSubgrid() == FALSE)
 
    /* Clean up this subgrid now that it is acceptable. */
 
    Subgrid->CleanUp();
 
    /* Go to the next grid in the queue. */
 
    index++;
 
  } // end: while (index < NumberOfSubgrids)
 
  return SUCCESS;
}
