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

  /* Return if not running Cluster SMBH feedback */
  if (ProblemType != 108)
    return SUCCESS;

  /* Return if we do not want to calculate the cold gas mass */
  if (ClusterSMBHCalculateGasMass == 0)
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
  MassUnits = DensityUnits*POW(LengthUnits,3);

  float ColdGasMassMsun=ClusterSMBHColdGasMass*MassUnits/SolarMass;
  if (ClusterSMBHCalculateGasMass >1 ){   //2-calculate&remove,3-calculate&remove&re-orient,4-Bondi
    if (ClusterSMBHCalculateGasMass == 3){
       ClusterSMBHJetDim = floor(Time*TimeUnits/(1.0e6*3.1557e7*ClusterSMBHJetPrecessionPeriod));   //ClusterSMBHJetPrecessionPeriod is now the Dim changing period.
    }
    if (ColdGasMassMsun  > 0.000001) {
       ClusterSMBHFeedbackSwitch = TRUE;
      if (ClusterSMBHCalculateGasMass == 4){  // no matter how much cold gas there is; accretiontime is now =dtFixed
         ClusterSMBHJetMdot = (ColdGasMassMsun/(ClusterSMBHAccretionTime))/2.0;  // AccretionTime already in s; Mdot in Msun/s. Devide it by 2 because Mdot is for only one jet.
         ClusterSMBHJetEdot = (ClusterSMBHAccretionEpsilon*ClusterSMBHJetMdot * SolarMass) * POW(clight,2)/1.0e44;   //for one jet
         } 
         else {
       ClusterSMBHJetMdot = (ColdGasMassMsun/(ClusterSMBHAccretionTime*1e6))/2.0;  // AccretionTime from Myr to yr; reset Mdot, still in Msun/yr. Devide it by 2 because Mdot is for only one jet.
       ClusterSMBHJetEdot = (ClusterSMBHAccretionEpsilon*ClusterSMBHJetMdot * SolarMass/3.1557e7) * POW(clight,2)/1.0e44;   //for one jet
       }
      }
    else
       ClusterSMBHFeedbackSwitch = FALSE;    // if there is not enough ColdGas, then do not turn jet on.
  } // end if ClusterSMBHCalculateGasMass > 1
  if (ClusterSMBHCalculateGasMass == 1) {
    int LastClusterSMBHFeedbackSwitch = ClusterSMBHFeedbackSwitch;
    if (ColdGasMassMsun < 1.0e5)
       ClusterSMBHFeedbackSwitch = FALSE;
    else
       ClusterSMBHFeedbackSwitch = (ColdGasMassMsun <= ClusterSMBHEnoughColdGas && LastClusterSMBHFeedbackSwitch == FALSE) ? FALSE : TRUE;

    if (LastClusterSMBHFeedbackSwitch == FALSE && ClusterSMBHFeedbackSwitch == TRUE) {
       ClusterSMBHStartTime = Time + ClusterSMBHTramp*1.0e6*3.1557e7/TimeUnits;
       if (ClusterSMBHJetPrecessionPeriod < 0.00001 & ClusterSMBHJetPrecessionPeriod > -0.00001)  //if precession off (set to 0), change the angle of the jets
       ClusterSMBHJetAnglePhi += 0.5;
    if (ClusterSMBHJetPrecessionPeriod < -0.00001)  //if precession negative, change the jet dimension
       ClusterSMBHJetDim += 1;
    }
  } // end if ClusterSMBHCalculateGasMass == 1
    

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *fptr=fopen("MT.out","a");
    fprintf(fptr,"Time, ClusterSMBHJetMdot, ClusterSMBHAccretionTime, and Total ClusterSMBHColdGasMass in Msun = %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n", Time, ClusterSMBHJetMdot, ClusterSMBHAccretionTime, ColdGasMassMsun);
    fclose(fptr);
  }
  return SUCCESS;
}
