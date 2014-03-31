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
#include "ErrorExceptions.h"
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
    ENZO_VFAIL("PE %"ISYM" NumberOfSubgrids > MAX_NUMBER_OF_SUBGRIDS in IdentifyNewSubgridsBySignature\n", MyProcessorNumber)
  }
 
  int index = 0;

  while (index < NumberOfSubgrids) {
 
    Subgrid = SubgridList[index];
 
    /* Shrink this subgrid (if necessary) to produce the smallest box. */
 
    if (Subgrid->ShrinkToMinimumSize() == FAIL) {
      ENZO_FAIL("Error in ProtoSubgrid->ShrinkToMinimumSize.");
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
	 ENZO_FAIL("Error in ProtoSubgrid->FindGridsByZeroSignature.");
	}
 
	/* Error check. */
 
	if (NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS) {
	  ENZO_FAIL("Increase MAX_NUMBER_OF_SUBGRIDS in IdentifyNewSubgridsBySignature.");
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

	    int GridRank = NewSubgrid->ReturnGridRank();
	    int DimLong = NewSubgrid->ReturnNthLongestDimension(0);
	    int DimShort = NewSubgrid->ReturnNthLongestDimension(GridRank-1);
	    
	    float ratio = float(DimLong)/float(DimShort);
	    if (ratio > 3.0)
	      printf("New subgrid has ratio = %g\n",ratio);

	  }

	  
	  break; // break out of the loop over dimensions
	}
	
      } // end: for (i = 0; i < MAX_DIMENSION; i++)
 
      /* If we couldn't partition by simple zeros then try the 2nd deriv. */
 
      if (NumberOfNewGrids <= 1) {

	int MaxZeroCrossingStrength = -1, StrongestDim = -1, TempInt;

	/* First check large axis ratio and split on that if necessary */
	
#define CRITICAL_RATIO 3.0


	dim = Subgrid->ReturnNthLongestDimension(0);
	Subgrid->LargeAxisRatioCheck(StrongestDim, GridEnds, CRITICAL_RATIO);

	int split_large_axis = 0;
	if (StrongestDim > -1){
	  
	  int *Dims;
	
	  split_large_axis = 1;
	  Dims = Subgrid->ReturnGridDimension();
	  //printf("StrongestDim = %d and Dims are %d, %d, %d\n",StrongestDim,Dims[0],Dims[1],Dims[2]);
	  //printf("New GridEnds: %d %d %d %d\n",GridEnds[StrongestDim*2][0],GridEnds[StrongestDim*2][1],GridEnds[StrongestDim*2+1][0],GridEnds[StrongestDim*2+1][1]);

	}
	if (StrongestDim == -1) {
 
	  /* Now Compute the zero crossings in the second derivaties of all the
	     signatures. */
 	
	  //	for (i = 0; i < MAX_DIMENSION; i++) {
	  for (i = 0; i < 1; i++) {
	
	    if ((dim = Subgrid->ReturnNthLongestDimension(i)) < 0)
	      break;
	
	    if (Subgrid->ComputeSecondDerivative(dim, TempInt,
						 &GridEnds[dim*2]) == FAIL) {
	      ENZO_FAIL("Error in ProtoSubgrid->ComputeSecondDerivative.\n");
	    }

	    int NewGridWidth_1, NewGridWidth_2;
	    NewGridWidth_1 = GridEnds[dim*2][1]-GridEnds[dim*2][0];
	    NewGridWidth_2 = GridEnds[dim*2+1][1]-GridEnds[dim*2+1][0];

	    if (TempInt > MaxZeroCrossingStrength and NewGridWidth_1 > MinimumSubgridEdge and NewGridWidth_2 > MinimumSubgridEdge) {
	      StrongestDim = dim;
	      MaxZeroCrossingStrength = TempInt;
	    }
	
	  } // end: for (i = 0; i < MAX_DIMENSION; i++)
	
	  /* Error check. */
	

	  if (StrongestDim < 0)
	    Subgrid->LargeAxisRatioCheck(StrongestDim, GridEnds, 0.0);

	  if (StrongestDim < 0) {
	    ENZO_FAIL("Error in IdentifyNewSubgridsBySignature.");
	  }

	} // end: if (StrongestDim == -1)


	//if (StrongestDim < 0)
	//  Subgrid->LargeAxisRatioCheck(StrongestDim, GridEnds, CRITICAL_RATIO);
	
	/* Create new subgrids (two). */

	int GridRank = Subgrid->ReturnGridRank();
	int DimLong = Subgrid->ReturnNthLongestDimension(0);
	int DimShort = Subgrid->ReturnNthLongestDimension(GridRank-1);
	
	int *Dims;
	Dims = Subgrid->ReturnGridDimension();

	if (StrongestDim != DimLong)
	  printf("StrongestDim is %d; DimLong is %d; DimShort is %d\n",StrongestDim, DimLong, DimShort);
	

	if (StrongestDim == DimLong){
	  int NewDimLong_1 = GridEnds[StrongestDim*2][1]-GridEnds[StrongestDim*2][0];
	  int NewDimLong_2 = GridEnds[StrongestDim*2+1][1]-GridEnds[StrongestDim*2+1][0];
	  
	  int DimMed =  Subgrid->ReturnNthLongestDimension(1);

	  int NewMaxDim_1 = max(NewDimLong_1,DimMed);
	  int NewMaxDim_2 = max(NewDimLong_2,DimMed);

	  float ratio_1 = float(NewMaxDim_1)/float(DimShort);
	  float ratio_2 = float(NewMaxDim_2)/float(DimShort);
	  
	  //printf("New grid ratios are %g and %g\n",ratio_1,ratio_2);
	}


	if (split_large_axis == 1)
	  printf("Splitting with New GridEnds: %d %d %d %d\n",GridEnds[StrongestDim*2][0],GridEnds[StrongestDim*2][1],GridEnds[StrongestDim*2+1][0],GridEnds[StrongestDim*2+1][1]);
	SubgridList[index] = new ProtoSubgrid;
	SubgridList[NumberOfSubgrids++] = new ProtoSubgrid;
	Subgrid->CopyToNewSubgrid(StrongestDim, GridEnds[StrongestDim*2][0],
				  GridEnds[StrongestDim*2][1],
				  SubgridList[index]);
	Subgrid->CopyToNewSubgrid(StrongestDim, GridEnds[StrongestDim*2+1][0],
				  GridEnds[StrongestDim*2+1][1],
				  SubgridList[NumberOfSubgrids-1]);

	
	int GridRank_1 = SubgridList[index]->ReturnGridRank();
        int DimLong_1 = SubgridList[index]->ReturnNthLongestDimension(0);
        int DimShort_1 = SubgridList[index]->ReturnNthLongestDimension(GridRank-1);

	int *Dims_1,*Dims_2;

	Dims_1 = SubgridList[index]->ReturnGridDimension();


	float ratio_1 = float(Dims_1[DimLong_1])/float(Dims_1[DimShort_1]);
	

	int GridRank_2 = SubgridList[NumberOfSubgrids-1]->ReturnGridRank();
        int DimLong_2 = SubgridList[NumberOfSubgrids-1]->ReturnNthLongestDimension(0);
        int DimShort_2 = SubgridList[NumberOfSubgrids-1]->ReturnNthLongestDimension(GridRank-1);
 
	Dims_2 = SubgridList[NumberOfSubgrids-1]->ReturnGridDimension();
	float ratio_2 = float(Dims_2[DimLong_2])/float(Dims_2[DimShort_2]);

	if (ratio_1 > 3.0 or ratio_2 > 3.0)
	  printf("New grid ratios are %g and %g; %d %d %d; %d %d %d; %d %d %d\n",ratio_1,ratio_2,Dims_1[0],Dims_1[1],Dims_1[2],Dims_2[0],Dims_2[1],Dims_2[2],Dims[0],Dims[1],Dims[2]);






	//if (split_large_axis == 1)
	printf("Breaking by zero-crossing. dim=%"ISYM"  break=%"ISYM"-%"ISYM"/%"ISYM"-%"ISYM"\n\n",
	       StrongestDim,
	       GridEnds[StrongestDim*2][0], GridEnds[StrongestDim*2][1],
	       GridEnds[StrongestDim*2+1][0], GridEnds[StrongestDim*2+1][1]);
	

      }
 
      /* Delete the old subgrid and set Subgrid to the (first) new grid. */
 
      delete Subgrid;
      Subgrid = SubgridList[index];
 
      /* Shrink this subgrid (if necessary) to produce the smallest box. */
 
      if (Subgrid->ShrinkToMinimumSize() == FAIL) {

	ENZO_FAIL("Error in ProtoSubgrid->ShrinkToMinimumSize.");
      }
 
    } // end: while (Subgrid->AcceptableSubgrid() == FALSE)
 
    /* Clean up this subgrid now that it is acceptable. */
 
    Subgrid->CleanUp();
 
    /* Go to the next grid in the queue. */
 
    index++;
 
  } // end: while (index < NumberOfSubgrids)
 
  return SUCCESS;
}
