/***********************************************************************
/
/  ADD ANY MASS MARKED FOR ACCRETION ONTO THE STAR PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);

int Star::Accrete(void)
{

  if (CurrentGrid == NULL)
    return SUCCESS;

  int dim, i, n, count;
  FLOAT time = CurrentGrid->Time;
  float dt = CurrentGrid->dtFixed;
  float this_dt, ratio1, ratio2, new_vel;
  
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, time);

  /* Sum up prescribed mass accretion, then add to mass */

  n = 0;
  while (n < naccretions && accretion_time[n] <= time+dt) {
    if (n+1 == naccretions)
      this_dt = time - accretion_time[n];
    else
      this_dt = accretion_time[n+1] - accretion_time[n];
    DeltaMass += accretion_rate[n++] * this_dt * TimeUnits;
  }
  Mass += DeltaMass;
  FinalMass += DeltaMass;

  /* Conserve momentum: change star particle velocity due to accreted
     material */

  ratio2 = DeltaMass / Mass;
  ratio1 = 1.0 - ratio2;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    vel[dim] = ratio1 * vel[dim] + ratio2 * delta_vel[dim];
    delta_vel[dim] = 0.0;
  }

  /* Remove these entries in the accretion table */

  count = 0;
  naccretions -= n;
  if (naccretions > 0) {
    FLOAT *temp_time = new FLOAT[naccretions];
    float *temp_rate = new float[naccretions];
    for (i = 0; i < naccretions+n; i++)
      if (accretion_time[i] <= time+dt) {
	temp_rate[count] = accretion_rate[i];
	temp_time[count] = accretion_time[i];
	count++;
      }
  
    delete [] accretion_rate;
    delete [] accretion_time;
    accretion_rate = temp_rate;
    accretion_time = temp_time;
  } else {  
    // No more accretion data
    accretion_rate = NULL;
    accretion_time = NULL;
  }
  
  return SUCCESS;
}
