/*------------------------------------------------------------------------
  SET EVOLVING STELLAR MASS THRESHOLD
  By Brian O'Shea

  History:
     23 April 2019 : BWO -- Created

  Note:
     StarMakerMinimumMassRamp = 0 is off
     StarMakerMinimumMassRamp = 1 is linear evolution of mass in time
     StarMakerMinimumMassRamp = 2 is linear evolution of mass in redshift
     StarMakerMinimumMassRamp = 3 is exponential evolution of mass in time
     StarMakerMinimumMassRamp = 4 is exponential evolution of mass in redshift

  If StarMakerMinimumMassRamp > 0, all values of the ramp parameters (starting and
  ending times and masses) are set in ReadParameterFile.C, and tests are made there 
  to ensure that the user has set them.

------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"

void my_exit(int status);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int SetStellarMassThreshold(FLOAT time)
{

  int timestep, i;
  FLOAT a, dadt, redshift=0.0;
  float early_mass, late_mass, current_mass, float_time=0.0, float_redshift=0.0;

  /* Return if not used */
  if (StarMakerMinimumMassRamp == 0)
    return SUCCESS;

  /* Calculate redshift */
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    redshift = (1 + InitialRedshift)/a - 1;
  }

  /* recast redshift and time to be the same precision as everything else; only important if
     the two different floating-point precisions are different. */
  float_redshift = (float) redshift;
  float_time = (float) time;

  if(StarMakerMinimumMassRamp == 1 || StarMakerMinimumMassRamp == 3){  // interpolation in time

    /* Set early and late masses in linear or log */
    if(StarMakerMinimumMassRamp == 1){ // mass evolution linear in time
      early_mass = StarMakerMinimumMassRampStartMass;
      late_mass = StarMakerMinimumMassRampEndMass;
    } else { // mass evolution exponential in time
      early_mass = log10(StarMakerMinimumMassRampStartMass);
      late_mass = log10(StarMakerMinimumMassRampEndMass);
    }

    /* set current stellar minimum mass threshold */
    if(float_time <= StarMakerMinimumMassRampStartTime){ // if time is before ramp start time, use early mass
      current_mass = early_mass;
    } else if (float_time >= StarMakerMinimumMassRampEndTime){ // if time is after ramp end time, use late mass
      current_mass = late_mass;
    } else {  // otherwise, linearly interpolate between start and end
      current_mass = early_mass + (float_time - StarMakerMinimumMassRampStartTime)
	* (late_mass-early_mass)/(StarMakerMinimumMassRampEndTime-StarMakerMinimumMassRampStartTime);  
    }

    /* set StarMakerMinimumMass correctly */
    if(StarMakerMinimumMassRamp == 1){
      StarMakerMinimumMass = current_mass;
    } else {
      StarMakerMinimumMass = POW(10.0,current_mass);
    }

  } else if(StarMakerMinimumMassRamp == 2 || StarMakerMinimumMassRamp == 4){  // interpolation in redshift

    /* set early and late masses in linear or log */
    if(StarMakerMinimumMassRamp == 2){ // mass evolution linear in redshift
      early_mass = StarMakerMinimumMassRampStartMass;
      late_mass = StarMakerMinimumMassRampEndMass;
    } else { // mass evolution exponential in time
      early_mass = log10(StarMakerMinimumMassRampStartMass);
      late_mass = log10(StarMakerMinimumMassRampEndMass);
    }

    /* set current stellar minimum mass threshold */
    if(float_redshift >= StarMakerMinimumMassRampStartTime){ // if redshift is prior to ramp start redshift, use early mass
      current_mass = early_mass;
    } else if (float_redshift <= StarMakerMinimumMassRampEndTime){ // if redshift is after ramp end redshift, use late mass
      current_mass = late_mass;
    } else {  // otherwise, linearly interpolate between start and end
      current_mass = early_mass + (float_redshift - StarMakerMinimumMassRampStartTime)
	* (late_mass-early_mass)/(StarMakerMinimumMassRampEndTime-StarMakerMinimumMassRampStartTime);
    }

    /* set StarMakerMinimumMass correctly */
    if(StarMakerMinimumMassRamp == 2){
      StarMakerMinimumMass = current_mass;
    } else {
      StarMakerMinimumMass = POW(10.0,current_mass);
    }

  } else {  // user has made a poor choice
    fprintf(stderr,"SetStellarMassThreshold:  StarMakerMinimumMassRamp improperly set!\n");
    my_exit(EXIT_FAILURE);
  }
  
  if(debug){
    printf("SetStellarMassThreshold:  StarMakerMinimumMass set to %"FSYM" at time %"PSYM" (redshift %"PSYM")\n",
	   StarMakerMinimumMass, time, redshift);
  }
  
  return SUCCESS;

}
