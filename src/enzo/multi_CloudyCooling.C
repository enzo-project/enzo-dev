/***********************************************************************
/
/  Solve Cloudy + MultiSpecies cooling.
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

int multi_CloudyCooling(float *density,float *totalenergy,float *gasenergy,
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
			float HeIShieldFactor,float HeIIShieldFactor)
{

  /* Compute size (in floats) of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Declare some variables.

  int i,iter,dtIterMax,cellsCompleted;
  float energy;
  float dtHalf,dtTake,dtTolerance;

  // Grid arrays.

  bool *iterationMask;
  float *edot, *timeElapsed;

  /* Set max iterations and tolerance. */

  dtIterMax = 50000;
  dtTolerance = 1e-10;
  dtHalf = dtFixed/2.;

  // Set up edot and iteration mask arrays.

  iterationMask = new bool[size];
  edot = new float[size];
  timeElapsed = new float[size];

  // Loop over all grid cells.

  for (i = 0;i < size;i++) {

    // Set all iteration masks to true.
    iterationMask[i] = 1;

    // Set time elapsed to zero.
    timeElapsed[i] = 0.0;

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

  // Call cooling solver until all cells have been integrated over a full dt.

  iter = 0;
  cellsCompleted = 0;

  while ((cellsCompleted < size) && (iter < dtIterMax)) {

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

      if (iterationMask[i]) {

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

	// Calculate timestep as minimum of 10% cooling time, half grid time-step, and time remaining.

	dtTake = min(fabs(0.1*energy/edot[i]),min(dtFixed-timeElapsed[i],dtHalf));

	// update total energy and gas energy

	totalenergy[i] += edot[i]*dtTake;

	if (DualEnergyFormalism) {
	  gasenergy[i] += edot[i]*dtTake;
	}

	// update time elapsed

	timeElapsed[i] += dtTake;

	// Determine if cell has completed integration.
	if (fabs(dtFixed-timeElapsed[i]) <= dtFixed*dtTolerance) {
	  iterationMask[i] = 0;
	  cellsCompleted++;
	}

      } // if (iterationMask[i])

    } // for (i = 0;i < size;i++)

    iter++;

  } // while ((cellsCompleted < size) && (iter < dtIterMax))

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

  if (iter >= dtIterMax) {
    fprintf(stderr,"Max iterations reached in multi_CloudyCooling.\n");
    return FAIL;
  }

  // Delete arrays.
  delete [] iterationMask;
  delete [] edot;
  delete [] timeElapsed;

  return SUCCESS;

}
