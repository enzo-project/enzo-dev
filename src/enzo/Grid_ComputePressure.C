/***********************************************************************
/
/  GRID CLASS (COMPUTE THE PRESSURE FIELD AT THE GIVEN TIME)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state.
 
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
#include "hydro_rk/EOS.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::ComputePressure(FLOAT time, float *pressure,
                          float MinimumSupportEnergyCoefficient)
{
 
  /* declarations */
 
  float density, gas_energy, total_energy, min_pressure;
  float velocity1, velocity2 = 0, velocity3 = 0;
  int i, size = 1;

  /* Error Check */
 
  if (time < OldTime || time > Time) {
    ENZO_FAIL("requested time is outside available range.\n");
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
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* If using Zeus_Hydro, then TotalEnergy is really GasEnergy so don't
     subtract the kinetic energy term. */

 
  float OneHalf = 0.5;
  if (HydroMethod == Zeus_Hydro)
    OneHalf = 0.0;
 
  /* Loop over the grid, compute the thermal energy, then the pressure,
     the timestep and finally the implied timestep. */
 
  /* special loop for no interpolate. */
 
  if (time == Time)
 
    for (i = 0; i < size; i++) {
 
      if( EquationOfState == 1 ){
	pressure[i] = IsothermalSoundSpeed *IsothermalSoundSpeed * BaryonField[DensNum][i];
      }else{
	if( EquationOfState == 0 )
	  total_energy  = BaryonField[TENum][i];
	density       = BaryonField[DensNum][i];
	velocity1     = BaryonField[Vel1Num][i];
	if (GridRank > 1 || HydroMethod == MHD_Li)
	  velocity2   = BaryonField[Vel2Num][i];
	if (GridRank > 2 || HydroMethod == MHD_Li)
	  velocity3   = BaryonField[Vel3Num][i];
	
	if (EOSType > 0){

	  /* If using polytropic EOS, calculate pressure directly from density */

	  float e, h, cs, dpdrho, dpde;
	  EOS(pressure[i], density, e, h, cs, dpdrho, dpde, EOSType, 0);
	  
	  /* also reset energy */
	  float kineticE = OneHalf * (
	    velocity1*velocity1 + velocity2*velocity2 + velocity3*velocity3
	  );
	  BaryonField[TENum][i] = e + kineticE;

	  if (HydroMethod == MHD_RK || UseMHDCT) {
	    float B2 = pow(BaryonField[B1Num][i],2)
	      + pow(BaryonField[B2Num][i],2)
	      + pow(BaryonField[B3Num][i],2);
	    BaryonField[TENum][i] += OneHalf * B2 / density;
	  }
	
	} else { 
	  if (DualEnergyFormalism == 0){ 
	    /* gas energy = E - 1/2 v^2. */
	    gas_energy    = total_energy - OneHalf*(velocity1*velocity1 +
						    velocity2*velocity2 +
						    velocity3*velocity3);

	    if (HydroMethod == MHD_RK || UseMHDCT) {
	      float B2 = pow(BaryonField[B1Num][i],2)
	        + pow(BaryonField[B2Num][i],2)
	        + pow(BaryonField[B3Num][i],2);
	      gas_energy -= OneHalf * B2 / density;
	    }
	  } else {
	    gas_energy = BaryonField[GENum][i];
	  }

	  pressure[i] = (Gamma - 1.0)*density*gas_energy;

	  if (pressure[i] < tiny_number)
	    pressure[i] = tiny_number;

      min_pressure =
        MinimumSupportEnergyCoefficient * (Gamma - 1.0) * density * density;

         if (pressure[i] < min_pressure)
          pressure[i] = min_pressure;

        }

      }//EquationOfState == 0

    } // end of loop
 
  else
 
    /* general case: */
 
    for (i = 0; i < size; i++) {
 
      if( EquationOfState == 1 ){
	pressure[i] = coef   *   BaryonField[DensNum][i] +
	  coefold*OldBaryonField[DensNum][i];
	pressure[i] *= IsothermalSoundSpeed*IsothermalSoundSpeed;

      }else{

	total_energy  = coef   *   BaryonField[TENum][i] +
	                coefold*OldBaryonField[TENum][i];
	density       = coef   *   BaryonField[DensNum][i] +
                        coefold*OldBaryonField[DensNum][i];
	velocity1     = coef   *   BaryonField[Vel1Num][i] +
                        coefold*OldBaryonField[Vel1Num][i];
	if (GridRank > 1 || UseMHDCT)
	  velocity2   = coef   *   BaryonField[Vel2Num][i] +
	                coefold*OldBaryonField[Vel2Num][i];
	if (GridRank > 2 || UseMHDCT)
	  velocity3   = coef   *   BaryonField[Vel3Num][i] +
	                coefold*OldBaryonField[Vel3Num][i];
 
	if (EOSType > 0) {
	
	  /* If using polytropic EOS, calculate pressure directly from density */
	  float e, h, cs, dpdrho, dpde;
	  EOS(pressure[i], density, e, h, cs, dpdrho, dpde, EOSType, 0);

	  /* also reset energy */
	  float kineticE = OneHalf * (
	    velocity1*velocity1 + velocity2*velocity2 + velocity3*velocity3
	  );
	  BaryonField[TENum][i] = e + kineticE;

	  if (HydroMethod == MHD_RK || UseMHDCT) {
	    float B2 = pow(BaryonField[B1Num][i],2) 
	      + pow(BaryonField[B2Num][i],2) + pow(BaryonField[B3Num][i],2);
	    BaryonField[TENum][i] += OneHalf * B2 / density;
	  }
	  
	} else {
	  /* gas energy = E - 1/2 v^2. */
	  if (DualEnergyFormalism == 0) {
	    float kineticE = OneHalf * (
	      velocity1*velocity1 + velocity2*velocity2 + velocity3*velocity3
	    );
	    gas_energy = total_energy - kineticE;

	    if (HydroMethod == MHD_RK || UseMHDCT) {
	      float B2 = pow(BaryonField[B1Num][i],2) 
	        + pow(BaryonField[B2Num][i],2) + pow(BaryonField[B3Num][i],2);
	      gas_energy -= OneHalf * B2 / density;
	    }
	  } else {
	    gas_energy =  coef * BaryonField[GENum][i]
	      + coefold*OldBaryonField[GENum][i];
	  }

	
	  pressure[i] = (Gamma - 1.0)*density*gas_energy;
	
	  if (pressure[i] < tiny_number)
	    pressure[i] = tiny_number;

      min_pressure =
        MinimumSupportEnergyCoefficient * (Gamma - 1.0) * density * density;

      if (pressure[i] < min_pressure)
        pressure[i] = min_pressure;


       }
      } //EquationOfState == 0

    } /* end of loop over cells */
 
  /* Correct for Gamma from H2. */
 
  if (MultiSpecies > 1) {
 
    float TemperatureUnits=1, number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(Gamma-1.0), x, Gamma1, temp;
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1; 
 
    /* Find Multi-species fields. */
 
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
      H2IINum, DINum, DIINum, HDINum;
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
 
    /* Find the temperature units. */
 
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.\n");
    }
 
    for (i = 0; i < size; i++) {
 
      number_density =
	  0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
		BaryonField[HeIIINum][i]                        ) +
	        BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
                BaryonField[DeNum][i];
 
      nH2 = 0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = tiny_number;
      temp = max(TemperatureUnits*pressure[i]/(number_density + nH2), 1);
 
      /* Only do full computation if there is a reasonable amount of H2.
	 The second term in GammaH2Inverse accounts for the vibrational
	 degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2/number_density > 1e-3) {
	x = 6100.0/temp;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      Gamma1 = 1.0 + (nH2 + number_density) /
	             (nH2*GammaH2Inverse + number_density * GammaInverse);
	
      /* Correct pressure with improved Gamma. */
 
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (MultiSpecies > 1)
 
  /* To emulate the opacity limit in turbulent star formation 
     simulations */
  
  float Gamma1 = Gamma;
  if ((ProblemType == 60 || ProblemType == 61) && SelfGravity == 1)

    for (i=0; i<size; i++) {
      Gamma1 = min(Gamma + (log10(BaryonField[DensNum][i])-8.0)*0.3999/2.5, 1.4);
      pressure[i] *= (Gamma1 - 1.0)/(Gamma - 1.0);
    }


  return SUCCESS;
}
