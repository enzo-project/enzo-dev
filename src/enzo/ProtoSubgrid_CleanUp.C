/***********************************************************************
/
/  PROTOSUBGRID CLASS (FREE UNNECESSARY MEMORY ALLOCATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1995
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
 
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
int ProtoSubgrid::CleanUp()
{
 
  delete [] GridFlaggingField;
  GridFlaggingField = NULL;
 
  /* Delete signatures unless dim=1 in which case Signature[0] = GFF */
 
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] Signature[dim];
    Signature[dim] = NULL;
  }
 
  /* Set the dimension according to the refinement and ghost zones. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    GridDimension[dim] = GridDimension[dim]*RefineBy + 2*NumberOfGhostZones;
  }
 
  /* Compute efficiency & report. */
 
  int i;
  if (debug1) {
    printf("ProtoSubgrid: efficiency = %6.1"FSYM"%% (%"ISYM"/%"ISYM") dims=",
	   float(NumberFlagged)/float(size)*100.0, NumberFlagged, size);
    for (i = 0; i < GridRank; i++)
      printf("%"ISYM" ", (GridDimension[i] - 2*NumberOfGhostZones)/RefineBy);
    printf("\n");
  }
 
 
  return SUCCESS;
}
