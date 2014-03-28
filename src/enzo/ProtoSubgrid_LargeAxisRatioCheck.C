/***********************************************************************
/
/  PROTOSUBGRID CLASS (CHECK FOR LARGE AXIS RATIO AND SPLIT IF TOO LARGE)
/
/  written by: Greg Bryan
/  date:       Mar, 2014
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int ProtoSubgrid::LargeAxisRatioCheck(int &SplitDim, int GridEnds[MAX_DIMENSION*2][2], float CriticalRatio)
{

  SplitDim = -1;

  if (GridRank == 1)
    return SUCCESS;

  /* Find longest and shortest dimension. */

  int DimLong = this->ReturnNthLongestDimension(0);
  int DimShort = this->ReturnNthLongestDimension(GridRank-1);

  /* If ratio of longest to shortest is above a critical value, then split grid
     in half along longest dimension. */

  if (float(GridDimension[DimLong])/float(GridDimension[DimShort]) > CriticalRatio) {
    SplitDim = DimLong;
    int Center = (GridDimension[SplitDim]-1)/2;
    GridEnds[SplitDim*2][0] = StartIndex[SplitDim];
    GridEnds[SplitDim*2][1] = StartIndex[SplitDim] + Center;
    GridEnds[SplitDim*2+1][0] = min(GridEnds[SplitDim*2][1]+1, EndIndex[SplitDim]);
    GridEnds[SplitDim*2+1][1] = EndIndex[SplitDim];
  }
 
  return SUCCESS;
}
