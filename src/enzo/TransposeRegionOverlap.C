/***********************************************************************
/
/  DETERMINE WHETHER TWO TRANSPOSE REGIONS OVERLAP
/
/  written by: Greg Bryan (Originally in CommunicationTranspose)
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"


int TransposeRegionOverlap(region *FromRegion, region *ToRegion, int i, int j,
			   region *Sends, region *Receives, 
			   int &nrecv, int &nsend,
			   int &SendSize, int &ReceiveSize)

{

  int dim, size;
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];

  /* Determine if there is an overlap. */

  for (dim = 0, size = 1; dim < MAX_DIMENSION; dim++) {
    LeftIndex[dim] = max(FromRegion[j].StartIndex[dim],
			 ToRegion[i].StartIndex[dim]);
    RightIndex[dim] = 
      min(FromRegion[j].StartIndex[dim] + FromRegion[j].RegionDim[dim],
	  ToRegion[i].StartIndex[dim] + ToRegion[i].RegionDim[dim])-1;
    size *= max(RightIndex[dim] - LeftIndex[dim] + 1, 0);
  } // ENDFOR dim

  /* If there is an overlap, add it to the list of sends/receives. */

  if (size > 0) {

    if (MyProcessorNumber == FromRegion[j].Processor) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	Sends[nsend].StartIndex[dim] = LeftIndex[dim] -
	  FromRegion[j].StartIndex[dim];
	Sends[nsend].RegionDim[dim] = RightIndex[dim] -
	  LeftIndex[dim] + 1;
      }
      SendSize += size;
      Sends[nsend++].Processor = j;
    }

    if (MyProcessorNumber == ToRegion[i].Processor) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	Receives[nrecv].StartIndex[dim] = LeftIndex[dim] -
	  ToRegion[i].StartIndex[dim];
	Receives[nrecv].RegionDim[dim] = RightIndex[dim] -
	  LeftIndex[dim] + 1;
      }
      ReceiveSize += size;
      Receives[nrecv++].Processor = i;
    }

  } // ENDIF size > 0

  return SUCCESS;

}
