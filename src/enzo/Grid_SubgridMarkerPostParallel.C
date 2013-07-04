/***********************************************************************
/
/  GRID CLASS (CONVERT STORED GRID IDS AND LEVELS TO SUBGRID MARKERS)
/
/  written by: John Wise
/  date:       May, 2010
/  modified1:  
/
/  PURPOSE:
/
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
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
#include "communication.h"

#define LEVEL_BIT_OFFSET 24
#define PACKED_INT_SIZE 32

/* function prototypes */

int grid::SubgridMarkerPostParallel(HierarchyEntry **Grids[], int *NumberOfGrids)
{

  /* Only process if this is the local processor, and the data wasn't
     originally on the processor before load balancing. */

  if (MyProcessorNumber != ProcessorNumber ||
      ProcessorNumber == OriginalProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  const int mask = (1 << LEVEL_BIT_OFFSET) - 1;
  int i, j, k, dim, index, size, GridID, GridLevel;
  long packed_int;
  
  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* check if the field has been allocated.  It should have been
     earlier! */ 

  if (SubgridMarker == NULL)
    ENZO_FAIL("SubgridMarker not allocated!");

  /* Convert packed integer (grid ID in bits 0-23 and level in bits
     24-31) into a grid pointer.  Do not overwrite pre-existing
     pointer because there will be overlap on the slabs. When
     addresses are 64-bit, we have to crop the address. */

  for (index = 0; index < size; index++) {
    if ((long) SubgridMarker[index] == INT_UNDEFINED) {
      SubgridMarker[index] = NULL;
    } else {
      packed_int = ((long) SubgridMarker[index] & 0xffffffff);
      GridID = packed_int & mask;
      GridLevel = packed_int >> LEVEL_BIT_OFFSET;
      if (GridID < 0 || GridID >= NumberOfGrids[GridLevel])
	ENZO_VFAIL("BBP%d: packed_int[%d] = %d, Grid %d, MarkerGrid/Level %d/%d\n", 
		   MyProcessorNumber, index, packed_int, this->ID, GridLevel, GridID);
      SubgridMarker[index] = Grids[GridLevel][GridID]->GridData;
    }
  } // ENDFOR index

  return SUCCESS;
  
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

int grid::SubgridMarkerPostParallelGZ(grid *Parent, HierarchyEntry **Grids[],
				      int *NumberOfGrids)
{

  /* Return if this grid is not on this processor or the parent is on
     the same processor. */

  if (MyProcessorNumber != ProcessorNumber ||
      ProcessorNumber == Parent->ProcessorNumber)
    return SUCCESS;

  /* declarations */
    
  const int mask = (1 << LEVEL_BIT_OFFSET) - 1;
  int i, j, k, dim, index, size, GridID, GridLevel;
  long packed_int;
  int GZStart[MAX_DIMENSION], GZEnd[MAX_DIMENSION];
  
  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  char *cellmask = new char[size];
  for (i = 0; i < size; i++) cellmask[i] = 0;

  /* check if the field has been allocated.  It should have been
     earlier! */ 

  if (SubgridMarker == NULL)
    ENZO_FAIL("SubgridMarker not allocated!");

  /* Convert packed integer (grid ID in bits 0-23 and level in bits
     24-31) into a grid pointer.  Do not overwrite pre-existing
     pointer because there will be overlap on the slabs. */

  for (dim = 0; dim < 3; dim++) {

    // Left face
    for (i = 0; i < 3; i++) {
      GZStart[i] = 0;
      GZEnd[i] = (i == dim) ? GridStartIndex[i]-1 : GridDimension[i]-1;
    }

    for (k = GZStart[2]; k <= GZEnd[2]; k++)
      for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	index = GRIDINDEX_NOGHOST(GZStart[0], j, k);
	for (i = GZStart[0]; i <= GZEnd[0]; i++, index++) {
	  if (cellmask[index] == 0) {
	    if ((long) SubgridMarker[index] == INT_UNDEFINED) {
	      SubgridMarker[index] = NULL;
	    } else {
	      // When addresses are 64-bit, we have to crop the address
	      packed_int = ((long) SubgridMarker[index] & 0xffffffff);
	      GridID = packed_int & mask;
	      GridLevel = packed_int >> LEVEL_BIT_OFFSET;
	      SubgridMarker[index] = Grids[GridLevel][GridID]->GridData;
	      cellmask[index] = 1;
	    }
	  } // ENDIF
	} // ENDFOR i
      } // ENDFOR j


    // Right face
    for (i = 0; i < 3; i++) {
      GZStart[i] = (i == dim) ? GridEndIndex[i]+1 : 0;
      GZEnd[i] = GridDimension[i]-1;
    }

    for (k = GZStart[2]; k <= GZEnd[2]; k++)
      for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	index = GRIDINDEX_NOGHOST(GZStart[0], j, k);
	for (i = GZStart[0]; i <= GZEnd[0]; i++, index++) {
	  if (cellmask[index] == 0) {
	    if ((long) SubgridMarker[index] == INT_UNDEFINED) {
	      SubgridMarker[index] = NULL;
	    } else {
	      // When addresses are 64-bit, we have to crop the address
	      packed_int = ((long) SubgridMarker[index] & 0xffffffff);
	      GridID = packed_int & mask;
	      GridLevel = packed_int >> LEVEL_BIT_OFFSET;
	      SubgridMarker[index] = Grids[GridLevel][GridID]->GridData;
	      cellmask[index] = 1;
	    }
	  } // ENDIF
	} // ENDFOR i
      } // ENDFOR j

  } // ENDFOR dim

  delete [] cellmask;
    
  return SUCCESS;
  
}
