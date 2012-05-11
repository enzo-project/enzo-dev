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

//int GetUnits(float *DensityUnits, float *LengthUnits,
//             float *TemperatureUnits, float *TimeUnits,
//             float *VelocityUnits, FLOAT Time);

int ClusterSMBHSumGasMass(HierarchyEntry *Grids[], int NumberOfGrids, int level)
{

  /* Return if we do not want to calculate the cold gas mass */
  if (ClusterSMBHCalculateGasMass != TRUE)
    return SUCCESS;

  int grid;

  /* Sum over all grids on this processor. */

  ClusterSMBHColdGasMass = 0;
  for (grid = 0; grid < NumberOfGrids; grid++) {
    Grids[grid]->GridData->ClusterSMBHEachGridGasMass(level);
  }

  /* Sum over all processors. */

  CommunicationAllSumValues(&ClusterSMBHColdGasMass, 1);

//  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
//    TimeUnits = 1.0, VelocityUnits = 1.0;
//  double MassUnits = 1.0;

//  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
//               &TimeUnits, &VelocityUnits, Time) == FAIL) {
//    fprintf(stderr, "Error in GetUnits.\n");
//    return FAIL;
//  }
//  MassUnits = DensityUnits*pow(LengthUnits,3);


  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Total ClusterSMBGColdGasMass = %g\n", ClusterSMBHColdGasMass); //*MassUnits/SolarMass);
  }

  return SUCCESS;
}
