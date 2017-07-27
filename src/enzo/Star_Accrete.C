/***********************************************************************
/
/  ADD ANY MASS MARKED FOR ACCRETION ONTO THE STAR PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             September, 2009
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::Accrete(void)
{
  if (UseSupernovaSeedFieldSourceTerms == 1){
    if (this->CurrentGrid == NULL ||(this->naccretions == 0))
      return SUCCESS;
  }
  else {
      if (this->CurrentGrid == NULL || 
          (this->naccretions == 0 && fabs(this->DeltaMass) < tiny_number))
      return SUCCESS;
  }

  const double Msun = 1.989e33, yr = 3.1557e7;
  int dim, i, n, count;
  FLOAT time = CurrentGrid->Time;
  float dt = CurrentGrid->dtFixed;
  float this_dt;
  
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  /* Sum up prescribed mass accretion, then add to mass */

  n = 0;
  while (n < naccretions && accretion_time[n] <= time+dt) {
    if (n+1 == naccretions)
      this_dt = time - accretion_time[n];
    else
      this_dt = accretion_time[n+1] - accretion_time[n];
    DeltaMass += accretion_rate[n++] * this_dt * TimeUnits;
  }
  
//  printf("star::Accrete: old_Mass = %lf, DeltaMass = %f\n", Mass, DeltaMass); 
  double old_mass = Mass;
  Mass += (double)(DeltaMass);
//  FinalMass += (double)(DeltaMass);
//  printf("star::Accrete: new_Mass = %lf, DeltaMass = %f\n", Mass, DeltaMass); 


  /* Conserve momentum: change star particle velocity due to accreted material */

  /* [1] For BlackHole, 
     it is now accurately done in Star_SubtractAccretedMassFromCell */

  /* [2] For Star Formation, 
     We can still do this in approximate way (JHW, Jan10) */

  double ratio1, ratio2, new_vel;

  if (type != MBH || type != BlackHole) {
    ratio2 = DeltaMass / Mass;
    ratio1 = 1.0 - ratio2;
    Metallicity = ratio1 * Metallicity + ratio2 * deltaZ;
    deltaZ = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      vel[dim] = ratio1 * vel[dim] + ratio2 * delta_vel[dim];
      delta_vel[dim] = 0.0;
    }
  } else {
    deltaZ = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = 0.0;
    }
  }

  /* [3] For MBH, 
     because gas mass is added to MBH from many cells with zero net momentum,
     just decrease the particle's velocity accordingly. */

  if (type == MBH) {
//    printf("star::Accrete: old_vel[1] = %g\n", vel[1]);
    vel[0] *= old_mass / Mass; 
    vel[1] *= old_mass / Mass;
    vel[2] *= old_mass / Mass; 
//    printf("star::Accrete: old_mass = %lf  ->  Mass = %lf\n", old_mass, Mass); 
//    printf("star::Accrete: new_vel[1] = %g\n", vel[1]);
  }


  /* Keep the last accretion_rate for future use */

  if (n > 0)  last_accretion_rate = accretion_rate[n-1]; 

  fprintf(stdout, "star::Accrete:  last_accretion_rate = %"GOUTSYM
	  " Msun/yr, time = %"GOUTSYM", "
	  "accretion_time[0] = %"GOUTSYM", this_dt = %"GOUTSYM
	  ", DeltaMass = %"GOUTSYM", Mass = %lf\n",
	  last_accretion_rate*yr, time, accretion_time[0], this_dt, DeltaMass, Mass);

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
