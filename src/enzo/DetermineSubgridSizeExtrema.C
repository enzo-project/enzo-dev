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

#define MINIMUM_EDGE 4
#define MINIMUM_SIZE 2000

int DetermineSubgridSizeExtrema(long_int NumberOfCells, int level, int MaximumStaticSubgridLevel)
{

  if (SubgridSizeAutoAdjust == FALSE)
    return SUCCESS;

  /* Now determine subgrid size parameters */

  int grids_per_proc = (level > MaximumStaticSubgridLevel) ?
    OptimalSubgridsPerProcessor : 8;

  MaximumSubgridSize = NumberOfCells / 
    (NumberOfProcessors * grids_per_proc);
  MaximumSubgridSize = max(MaximumSubgridSize, MINIMUM_SIZE);
  MinimumSubgridEdge = nint(pow(MaximumSubgridSize, 0.33333) * 0.25);
  MinimumSubgridEdge += MinimumSubgridEdge % 2;
  MinimumSubgridEdge = max(MinimumSubgridEdge, MINIMUM_EDGE);

  if (debug)
    printf("DetermineSGSize: MaxSubgridSize = %"ISYM", MinSubgridEdge = %"
	   ISYM", ncells = %lld\n",
	   MaximumSubgridSize, MinimumSubgridEdge, NumberOfCells);

  return SUCCESS;

}
