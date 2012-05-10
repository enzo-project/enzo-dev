/***********************************************************************
/
/  LOOPS OVER GRID AND SUMS COLD GAS COMPONENT
/
/  written by: Yuan Li
/  date:       May 2012
/  modified1: 
/
/  NOTES:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
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
#include "CommunicationUtilities.h"


extern float ClusterSMBHColdGasMass;  //yuan

int ClusterSMBHSumGasMass(HierarchyEntry *Grids[], int NumberOfGrids, int level)
{
   int grid;

   /* Sum over all grids on this processor. */

   ClusterSMBHColdGasMass = 0;
   for (grid = 0; grid < NumberOfGrids; grid++) {
      Grids[grid]->GridData->ClusterSMBHEachGridGasMass(level);
   }

   /* Sum over all processors. */

   CommunicationAllSumValues(&ClusterSMBHColdGasMass, 1);

   if (MyProcessorNumber == ROOT_PROCESSOR) {
     printf("Total ClusterSMBGColdGasMass = %g\n", ClusterSMBHColdGasMass);
   }

   return SUCCESS;
}
