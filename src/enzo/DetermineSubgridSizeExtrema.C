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

int DetermineSubgridSizeExtrema(LevelHierarchyEntry *LevelArray[],
				int level)
{

  if (SubgridSizeAutoAdjust == FALSE)
    return SUCCESS;

  int dim, Rank, size;
  int NumberOfCells = 0;
  LevelHierarchyEntry *Temp;

  int Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  /* Sum number of cells on local grids */

  for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
    if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      Temp->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
      for (dim = 0, size = 1; dim < Rank; dim++)
	size *= Dims[dim];
      NumberOfCells += size;
    } // ENDIF local

  /* Sum them across processors */

  CommunicationAllSumValues(&NumberOfCells, 1);

  /* Now determine subgrid size parameters */

  MaximumSubgridSize = NumberOfCells / NumberOfProcessors / 
    OptimalSubgridsPerProcessor;
  MaximumSubgridSize = max(MaximumSubgridSize, MINIMUM_SIZE);
  MinimumSubgridEdge = nint(pow(MaximumSubgridSize, 0.33333) * 0.25);
  MinimumSubgridEdge += MinimumSubgridEdge % 2;
  MinimumSubgridEdge = max(MinimumSubgridEdge, MINIMUM_EDGE);

//  if (debug)
//    printf("DetermineSGSize: MaxSubgridSize = %"ISYM", MinSubgridEdge = %"
//	   ISYM", ncells = %"ISYM"\n",
//	   MaximumSubgridSize, MinimumSubgridEdge, NumberOfCells);

  return SUCCESS;

}
