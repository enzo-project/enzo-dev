/***********************************************************************
/
/  PROTOSUBGRID CLASS (DIVIDE THIS GRID BY FINDING ZEROS IN THE SIGNATURE)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int ProtoSubgrid::FindGridsByZeroSignature(int dim, int &NumberOfNewGrids,
				     int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2])
{
  /* Error check */
 
  if (dim >= GridRank) {
    fprintf(stderr, "Passed dim(%"ISYM") > GridRank(%"ISYM")\n", dim, GridRank);
    return FAIL;
  }
 
  if (Signature[dim] == NULL) {
    fprintf(stderr, "Signature %"ISYM" not yet computed.\n", dim);
    return FAIL;
  }
 
  /* Initialize */
 
  int i = 0;
  NumberOfNewGrids = 0;
 
  /* Loop over signature. */
 
  while (i < GridDimension[dim]) {
 
    /* Look for the start of a new subgrid. */
 
    if (Signature[dim][i] != 0) {
      GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + i;
 
      /* Now find the end of the subgrid. */
 
      while (i < GridDimension[dim] && Signature[dim][i] != 0)
	i++;
      GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + i-1;

      if ( NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS ) {
        fprintf(stderr, "PE %"ISYM" NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS in ProtoSubgrid_FindGridsByZeroSignature\n", MyProcessorNumber);
        return FAIL;
      }

    }
 
    /* Next zone in signature. */
 
    i++;
  }
 
  return SUCCESS;
}
