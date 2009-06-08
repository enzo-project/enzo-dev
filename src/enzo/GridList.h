/***********************************************************************
/
/  GRID LIST STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifndef GRID_LIST_DEFINED__
#define GRID_LIST_DEFINED__

#include <stdio.h>
#include "macros_and_parameters.h"

struct GridList
{
  GridList *NextGrid;
  int GridRank;
  int GridDimension[MAX_DIMENSION];
  FLOAT GridLeftEdge[MAX_DIMENSION];
  FLOAT GridRightEdge[MAX_DIMENSION];

  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];
  int NumberFlagged;
};


#endif
