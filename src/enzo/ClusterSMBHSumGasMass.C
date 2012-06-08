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
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

float ClusterSMBHColdGasMass;

int ClusterSMBHSumGasMass(HierarchyEntry *Grids[], int NumberOfGrids, int level)
{

  /* Return if we do not want to calculate the cold gas mass */
  if (ClusterSMBHCalculateGasMass != TRUE)
    return SUCCESS;

  /* Return if not on most-refined level. */
  if (level != MaximumRefinementLevel)
    return SUCCESS;

  int grid;

  /* Sum over all grids on this processor. */

  ClusterSMBHColdGasMass = 0;
  for (grid = 0; grid < NumberOfGrids; grid++) {
//    if (level == MaximumRefinementLevel)
      Grids[grid]->GridData->ClusterSMBHEachGridGasMass(level);
  }

  /* Sum over all processors. */

  CommunicationAllSumValues(&ClusterSMBHColdGasMass, 1);
  FLOAT Time = Grids[0]->GridData->ReturnTime();
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits = 1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  MassUnits = DensityUnits*pow(LengthUnits,3);

  float ColdGasMassMsun=ClusterSMBHColdGasMass*MassUnits/SolarMass;

  int LastClusterSMBHFeedbackSwitch = ClusterSMBHFeedbackSwitch;
  if (ColdGasMassMsun < 1.0e5)
    ClusterSMBHFeedbackSwitch = FALSE;
  else
    ClusterSMBHFeedbackSwitch = (ColdGasMassMsun <= ClusterSMBHEnoughColdGas && LastClusterSMBHFeedbackSwitch == FALSE) ? FALSE : TRUE;

  if (LastClusterSMBHFeedbackSwitch == FALSE && ClusterSMBHFeedbackSwitch == TRUE) {
    ClusterSMBHStartTime = Time + ClusterSMBHTramp*1.0e6*3.1557e7/TimeUnits;
    if (ClusterSMBHJetPrecessionPeriod < 0.00001)  //if precession off, change the angle of the jets
      ClusterSMBHJetAnglePhi += 0.5;
  }

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Time, ClusterSMBHStartTime, Switch, and Total ClusterSMBGColdGasMass in Msun = %g %g %d %g \n", Time, ClusterSMBHStartTime, ClusterSMBHFeedbackSwitch, ColdGasMassMsun);
  }

  return SUCCESS;
}
