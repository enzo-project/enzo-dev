/***********************************************************************
/
/  GRID CLASS (COMPUTE THE TEMPERATURE AT THE GIVEN TIME USING
/              GADGET EQUILIBRIUM COOLING)
/
/  written by: Brian O'Shea
/  date:       February 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
/  NOTE:  This routine is pretty blatantly borrowed from 
/         Grid_ComputePressure - if there are bugs here, check there as
/         well.
/
************************************************************************/

// compute the temperature here at the requested time, using
// Gadget's cooling code

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

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::GadgetComputeTemperature(FLOAT time, float *temperature)
{

  /* declarations */

  float density, gas_energy, total_energy;
  float velocity1, velocity2 = 0.0, velocity3 = 0.0;

  float *ne_guess,guess;

  // set electron fraction guess to be 1% for the first one, after that we just use
  // the previous cell
  guess=0.01;
  ne_guess = &guess;


  int i, size = 1;

  /* Error Check */

  if (time < OldTime || time > Time) {
    ENZO_FAIL("requested time is outside available range.");
  }

  /* Compute interpolation coefficients. */

  float coef, coefold;
  if (Time - OldTime > 0)
    coef    = (time - OldTime)/(Time - OldTime);
  else
    coef    = 1;

  coefold = 1 - coef;

  /* Compute the size of the grid. */

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  // get physical units
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1, 
    TemperatureUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  

  /* If using Zeus_Hydro, then TotalEnergy is really GasEnergy so don't
     subtract the kinetic energy term. */

  float OneHalf = 0.5;
  if (HydroMethod == Zeus_Hydro)
    OneHalf = 0.0;

  /* Loop over the grid, compute the thermal energy, then the temperature,
     the timestep and finally the implied timestep. */

  /* special loop for no interpolate. */

  if (time == Time)

    for (i = 0; i < size; i++) {

      total_energy  = BaryonField[TENum][i];
      density       = BaryonField[DensNum][i];
      velocity1     = BaryonField[Vel1Num][i];
      if (GridRank > 1)
	velocity2   = BaryonField[Vel2Num][i];
      if (GridRank > 2)
	velocity3   = BaryonField[Vel3Num][i];

      /* gas energy = E - 1/2 v^2. */
      gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
					      velocity2*velocity2 +
					      velocity3*velocity3);


      // convert energy, density to CGS
      gas_energy *= (VelocityUnits*VelocityUnits);

      density *= DensityUnits;

      temperature[i] = Gadgetconvert_u_to_temp(gas_energy, density, ne_guess);

      if (temperature[i] < 1.0)
	temperature[i] = 1.0;

    } // end of loop

  else

    /* general case: */

    for (i = 0; i < size; i++) {

      total_energy  = coef   *   BaryonField[TENum][i] + 
	              coefold*OldBaryonField[TENum][i];
      density       = coef   *   BaryonField[DensNum][i] + 
                      coefold*OldBaryonField[DensNum][i];
      velocity1     = coef   *   BaryonField[Vel1Num][i] + 
                      coefold*OldBaryonField[Vel1Num][i];

      if (GridRank > 1)
	velocity2   = coef   *   BaryonField[Vel2Num][i] + 
	              coefold*OldBaryonField[Vel2Num][i];
      if (GridRank > 2)
	velocity3   = coef   *   BaryonField[Vel3Num][i] + 
	              coefold*OldBaryonField[Vel3Num][i];


      /* gas energy = E - 1/2 v^2. */    
      gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
					      velocity2*velocity2 +
					      velocity3*velocity3);
      
      // convert energy, density to CGS
      gas_energy *= (VelocityUnits*VelocityUnits);

      density *= DensityUnits;

      temperature[i] = Gadgetconvert_u_to_temp(gas_energy, density, ne_guess);

      if (temperature[i] < 1.0)
	temperature[i] = 1.0;

    }


  return SUCCESS;
}
