/***********************************************************************
/
/  GRID CLASS (RETURN SUBGRID MARKER)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::ReturnSubgridMarker(int &cindex, FLOAT x, FLOAT y, FLOAT z)
{

  int gi, gj, gk;

  cindex = 0;

  if (x >= GridLeftEdge[0] && x <= GridRightEdge[0] &&
      y >= GridLeftEdge[1] && y <= GridRightEdge[1] &&
      z >= GridLeftEdge[2] && z <= GridRightEdge[2]) {

    gi = nint(x - GridLeftEdge[0]/CellWidth[0][0]) + GridStartIndex[0];
    gj = nint(y - GridLeftEdge[1]/CellWidth[1][0]) + GridStartIndex[1];
    gk = nint(z - GridLeftEdge[2]/CellWidth[2][0]) + GridStartIndex[2];

    cindex = GRIDINDEX(gi, gj, gk);
   
  }

  return SUCCESS;

}
