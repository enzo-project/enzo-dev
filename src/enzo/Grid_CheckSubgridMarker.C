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

/* function prototypes */

int grid::CheckSubgridMarker(void)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int i, j, k, index, dim, rank, GridDims[MAX_DIMENSION];
  FLOAT le[MAX_DIMENSION], re[MAX_DIMENSION], CellCenter[MAX_DIMENSION], DomainWidth[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  index = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, index++) {
	SubgridMarker[index]->ReturnGridInfo(&rank, GridDims, le, re);
	CellCenter[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];
	CellCenter[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];
	CellCenter[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k];
	// Not necessary to check the case where the cell center
	// equals the domain edge because the cell centers are always
	// offset by 1/2 of a cell width.
	for (dim = 0; dim < GridRank; dim++) {
	  if (CellCenter[dim] < DomainLeftEdge[dim])
	    CellCenter[dim] += DomainWidth[dim];
	  else if (CellCenter[dim] > DomainRightEdge[dim])
	    CellCenter[dim] -= DomainWidth[dim];
	}
	if (SubgridMarker[index]->PointInGridNB(CellCenter) == FALSE) {
	  printf("ThisLeftEdge  = %15.12f  %15.12f  %15.12f\n", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
	  printf("ThisRightEdge = %15.12f  %15.12f  %15.12f\n", GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
	  printf("SBLeftEdge    = %15.12f  %15.12f  %15.12f\n", le[0], le[1], le[2]);
	  printf("SBRightEdge   = %15.12f  %15.12f  %15.12f\n", re[0], re[1], re[2]);
	  printf("CellCenter    = %15.12f  %15.12f  %15.12f\n", CellCenter[0], CellCenter[1], CellCenter[2]);
	  ENZO_FAIL("");
	}
      }
    }
  }

  return SUCCESS;

}
