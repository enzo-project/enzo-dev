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
  
  If StarMakerMinimumMassRamp > 0, all values of the ramp parameters are set in
  ReadParameterFile.C and tests are made to ensure that they're not the standard
  values!

------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

void my_exit(int status);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int SetStellarMassThreshold(FLOAT time)
{

  int timestep, i;
  FLOAT a, dadt, redshift;
  float early_mass, late_mass, current_mass;

  /* Return if not used */
  if (StarMakerMinimumMassRamp == 0)
    return SUCCESS;

  /* If TimeType is redshift, calculate redshift */

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    redshift = (1 + InitialRedshift)/a - 1;
  }

  if(StarMakerMinimumMassRamp == 1 || StarMakerMinimumMassRamp == 3){  // interpolation in time

    /* set early and late masses in linear or log */
    if(StarMakerMinimumMassRamp == 1){ // mass evolution linear in time
      early_mass = StarMakerMinimumMassRampStartMass;
      late_mass = StarMakerMinimumMassRampEndMass;
    } else { // mass evolution exponential in time
      early_mass = log10(StarMakerMinimumMassRampStartMass);
      late_mass = log10(StarMakerMinimumMassRampEndMass);
    }

    /* set current stellar minimum mass threshold */
    if(time < StarMakerMinimumMassRampStartTime){
      current_mass = early_mass;
    } else if (time > StarMakerMinimumMassRampEndTime){
      current_mass = late_mass;
    } else {
      current_mass = early_mass + (time - StarMakerMinimumMassRampStartTime)
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
    if(redshift > StarMakerMinimumMassRampStartTime){
      current_mass = early_mass;
    } else if (redshift < StarMakerMinimumMassRampEndTime){
      current_mass = late_mass;
    } else {
      current_mass = early_mass + (redshift - StarMakerMinimumMassRampStartTime)
	* (late_mass-early_mass)/(StarMakerMinimumMassRampEndTime-StarMakerMinimumMassRampStartTime);
    }

    /* set StarMakerMinimumMass correctly */
    if(StarMakerMinimumMassRamp == 3){
      StarMakerMinimumMass = current_mass;
    } else {
      StarMakerMinimumMass = POW(10.0,current_mass);
    }

  } else {  // user has made a poor choice
    fprintf(stderr,"SetStellarMassThreshold:  StarMakerMinimumMassRamp improperly set!\n");
    my_exit(EXIT_FAILURE);
  }

  if(debut){
    printf("SetStellarMassThreshold:  StarMakerMinimumMass set to %"FSYM" at time %"PSYM"\n",
	   StarMakerMinimumMass,time);
  }
  
  return SUCCESS;

}
