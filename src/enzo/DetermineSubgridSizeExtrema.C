/***********************************************************************
/
/  DETERMINE MINIMUM AND MAXIMUM SIZES FOR SUBGRIDS
/
/  written by: John Wise
/  date:       August, 2010
/  modified1:
/
/  PURPOSE:  This is done for every level>0.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
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
#include "CommunicationUtilities.h"

// Minimum values to be used:

#define MINIMUM_EDGE 2
#define MINIMUM_SIZE 64

// Default values to be used, if no good estimate is availabke:

#define DEFAULT_MAXIMUM_SUBGRID_SIZE 2000
#define DEFAULT_MINIMUM_SUBGRID_EDGE 4

int DetermineSubgridSizeExtrema(long_int NumberOfCells, int level, int MaximumStaticSubgridLevel)
{

  if (SubgridSizeAutoAdjust == FALSE)
    return SUCCESS;

  /* Now determine subgrid size parameters */

  int grids_per_proc = (level > MaximumStaticSubgridLevel) ?
    OptimalSubgridsPerProcessor : 2;

  MaximumSubgridSize = NumberOfCells / 
    (NumberOfProcessors * grids_per_proc);
  MaximumSubgridSize = max(MaximumSubgridSize, MINIMUM_SIZE);
  MinimumSubgridEdge = nint(pow(MaximumSubgridSize, 0.33333) * 0.25);
  MinimumSubgridEdge += MinimumSubgridEdge % 2;
  MinimumSubgridEdge = max(MinimumSubgridEdge, MINIMUM_EDGE);

  /* If the NumberOfCells is zero, we have no good estimate to use to calculate
     the optimal values, so use some defaults, defined earlier. */

  if (NumberOfCells == 0) {
    MaximumSubgridSize = DEFAULT_MAXIMUM_SUBGRID_SIZE;
    MinimumSubgridEdge = DEFAULT_MINIMUM_SUBGRID_EDGE;
    printf("DetermineSGSize: Warning: NumberOfCells is 0, using defaults.\n");
  }

  if (debug)
    printf("DetermineSGSize: MaxSubgridSize = %"ISYM", MinSubgridEdge = %"
	   ISYM", ncells = %"ISYM"\n",
	   MaximumSubgridSize, MinimumSubgridEdge, NumberOfCells);

  return SUCCESS;

}
