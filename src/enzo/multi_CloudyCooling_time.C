/***********************************************************************
/
/  Solve Cloudy + MultiSpecies cooling time.
/
/  written by: Britton Smith
/  date:       March, 2006
/  modified1: 
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int SolveCloudyCooling(float *density,float *totalenergy,float *gasenergy,
		       float *velocity1,float *velocity2,float *velocity3,
		       float *De,float *HI,float *HII,
		       float *HeI,float *HeII,float *HeIII,
		       float *HM,float *H2I,float *H2II,
		       float *DI,float *DII,float *HDI,
		       float *metalDensity,
		       int *GridDimension,int GridRank,float dtFixed,
		       float afloat,float TemperatureUnits,float LengthUnits,
		       float aUnits,float DensityUnits,float TimeUnits,
		       int RadiationShield,float HIShieldFactor,
		       float HeIShieldFactor,float HeIIShieldFactor,
		       bool *iterationMask,float *edot);

int multi_CloudyCooling_time(float *density,float *totalenergy,float *gasenergy,
			     float *velocity1,float *velocity2,float *velocity3,
			     float *De,float *HI,float *HII,
			     float *HeI,float *HeII,float *HeIII,
			     float *HM,float *H2I,float *H2II,
			     float *DI,float *DII,float *HDI,
			     float *metalDensity,float *cooling_time,
			     int *GridDimension,int GridRank,float dtFixed,
			     float afloat,float TemperatureUnits,float LengthUnits,
			     float aUnits,float DensityUnits,float TimeUnits,
			     int RadiationShield,float HIShieldFactor,
			     float HeIShieldFactor,float HeIIShieldFactor)
{

  /* Compute size (in floats) of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Declare some variables.

  int i;
  float energy;
  float dtHalf,dtTake,dtTolerance;

  // Grid arrays.

  bool *iterationMask;
  float *edot;

  // Set up edot and iteration mask arrays.

  iterationMask = new bool[size];
  edot = new float[size];

  // Loop over all grid cells.

  for (i = 0;i < size;i++) {

    // Set all iteration masks to true.
    iterationMask[i] = 1;

    // convert density from comoving to proper

    density[i] /= pow(afloat,3);
    if (MultiSpecies > 0) {
      De[i]      /= pow(afloat,3);
      HI[i]      /= pow(afloat,3);
      HII[i]     /= pow(afloat,3);
      HeI[i]     /= pow(afloat,3);
      HeII[i]    /= pow(afloat,3);
      HeIII[i]   /= pow(afloat,3);
    }
    if (MultiSpecies > 1) {
      HM[i]      /= pow(afloat,3);
      H2I[i]     /= pow(afloat,3);
      H2II[i]    /= pow(afloat,3);
    }
    if (MultiSpecies > 2) {
      DI[i]      /= pow(afloat,3);
      DII[i]     /= pow(afloat,3);
      HDI[i]     /= pow(afloat,3);
    }
    if (CloudyCoolingData.CloudyCoolingGridRank >= 2) {
      metalDensity[i] /= pow(afloat,3);
    }

  }

  // Call cooling solver.

  if (SolveCloudyCooling(density,totalenergy,gasenergy,velocity1,velocity2,velocity3,
			 De,HI,HII,HeI,HeII,HeIII,HM,H2I,H2II,DI,DII,HDI,metalDensity,
			 GridDimension,GridRank,dtFixed,afloat,TemperatureUnits,LengthUnits,
			 aUnits,DensityUnits,TimeUnits,RadiationShield,HIShieldFactor,
			 HeIShieldFactor,HeIIShieldFactor,iterationMask,edot) == FAIL) {
    fprintf(stderr,"Error in SolveRadiativeCooling.\n");
    return FAIL;
  }

  for (i = 0;i < size;i++) {

    // Calculate energy.

    // Zeus - total energy is really gas energy
    if (HydroMethod == 2) {
      energy = totalenergy[i];
    }

    // PPM
    else {
      // with Dual Energy Formalism
      if (DualEnergyFormalism) {
	energy = gasenergy[i];
      }
      // without Dual Energy Formalism
      else {
	energy = totalenergy[i] - 0.5*velocity1[i]*velocity1[i];
	if (GridRank > 1)
	  energy -= 0.5*velocity2[i]*velocity2[i];
	if (GridRank > 2)
	  energy -= 0.5*velocity3[i]*velocity3[i];
      }
    }

    // Get cooling time.

    cooling_time[i] = -energy / edot[i];

  } // for (i = 0;i < size;i++)

  // convert density back to comoving

  for (i = 0;i < size;i++) {

    density[i] *= pow(afloat,3);
    if (MultiSpecies > 0) {
      De[i]      *= pow(afloat,3);
      HI[i]      *= pow(afloat,3);
      HII[i]     *= pow(afloat,3);
      HeI[i]     *= pow(afloat,3);
      HeII[i]    *= pow(afloat,3);
      HeIII[i]   *= pow(afloat,3);
    }
    if (MultiSpecies > 1) {
      HM[i]      *= pow(afloat,3);
      H2I[i]     *= pow(afloat,3);
      H2II[i]    *= pow(afloat,3);
    }
    if (MultiSpecies > 2) {
      DI[i]      *= pow(afloat,3);
      DII[i]     *= pow(afloat,3);
      HDI[i]     *= pow(afloat,3);
    }
    if (CloudyCoolingData.CloudyCoolingGridRank >= 2) {
      metalDensity[i] *= pow(afloat,3);
    }

  } // for (i = 0;i < size;i++)

  // Delete arrays.
  delete [] iterationMask;
  delete [] edot;

  return SUCCESS;

}
