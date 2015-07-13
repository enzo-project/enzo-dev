/***********************************************************************
/
/  GRID CLASS (Compute thermal conduction time-scale)
/
/  written by:  David A. Ventimiglia & Brian O'Shea
/  date:        December, 2009
/
/  PURPOSE:  Calculates the shortest time scale for thermal conduction
/  on a given grid patch.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"


// Function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

// Member functions
float grid::ComputeConductionTimeStep (float &dt) {

  dt = huge_number;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
  this->DebugCheck("ComputeConductionTimeStep");

  // Some locals
  int DensNum, TENum, GENum, Vel1Num, Vel2Num, Vel3Num;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0, aUnits = 1.0;
  FLOAT a = 1.0, dadt;
  double MassUnits = 1.0;
  float *rho;
  float dt_est, light_cross_time;
  double all_units;

  float SpitzerFraction;
  if (IsotropicConduction && AnisotropicConduction) {
    SpitzerFraction = max(IsotropicConductionSpitzerFraction, 
			  AnisotropicConductionSpitzerFraction);
  }
  else if (IsotropicConduction) {
    SpitzerFraction = IsotropicConductionSpitzerFraction;
  }
  else if (AnisotropicConduction) {
    SpitzerFraction = AnisotropicConductionSpitzerFraction;
  }
  else {
    return SUCCESS;
  }

  int size = 1, grid_index, right_side_index; 
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  float *Temp = new float[size];
  FLOAT dx = CellWidth[0][0];

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }

  // for conduction saturation
  double saturation_factor = 4.874e-20 / (DensityUnits * LengthUnits * dx);
                                        // 4.2 * lambda_e * mH
                                        // lambda_e from Jubelgas ea 2004
                                        // mH for converting rho into n_e
                                        // dx for dT/dx

  // get field identifiers
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  // get temperature
  if (this->ComputeTemperatureField(Temp) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  // mask for baryon density
  rho = BaryonField[DensNum];

  // Set up a struct to hold properties defined on cell faces
  struct cellface {float T, dT, kappa, dedt, rho;} l, r, cfzero;

  // zero struct
  cfzero.T = cfzero.dT = cfzero.kappa = cfzero.dedt = cfzero.rho = 0.0;

  // Find shortest time scale on the grid patch
  int GridStart[] = {0, 0, 0}, 
    GridEnd[] = {0, 0, 0};
  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;}

  dt = huge_number;

  /* timestep is calculated as dt < 0.5 * dx^2 / alpha, where
     alpha = thermal diffusivity = Kappa/(number density * k_boltz).
     The 0.5 is actually put in ComputeTimeStep 
     (as ConductionCourantSafetyFactor), and the rest is calculated
     here.  We first calculate the density and temperature parts in
     the loop (this is the only stuff that varies), and then multiply
     by all of the constants necessary to put this into Enzo internal
     units.  Since the saturation coefficient requires the temperature 
     gradients, this calculation now looks a lot like what is done in 
     Grid_ComputeHeat. */

  if (GridRank>0) {
    for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
      for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
	r = cfzero;
	for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i+1,j,k);

	  if(i == GridEnd[0]){
	    r = cfzero;
	  } else {

	    // get temperature, temperature gradient on + face of cell
	    // (the 'l' struct has it on the right face)
	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];

	    // kappa is the spitzer conductivity, which scales as 
	    // the temperature to the 2.5 power
	    r.kappa = POW(r.T, 2.5);
	    // conduction saturation
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	  }

	  dt_est = rho[ELT(i,j,k)] / r.kappa;
	  dt = min(dt, dt_est);
	}
      }
    }
  } // if (GridRank>0)

  if (GridRank>1) {
    for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
      for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
	r = cfzero;
	for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j+1,k);

	  if(j==GridEnd[1]){
	    r = cfzero;
	  } else {

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];
	    
	    r.kappa = POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));
	  }

	  dt_est = rho[ELT(i,j,k)] / r.kappa;
	  dt = min(dt, dt_est);
	}
      }
    }
  } // if (GridRank>1)
  
  if (GridRank>2) {
    for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
      for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
	r = cfzero;
	for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j,k+1);

	  if(k==GridEnd[2]){
	    r = cfzero;
	  } else {

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];

	    r.kappa = POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));
	  }

	  dt_est = rho[ELT(i,j,k)] / r.kappa;
	  dt = min(dt, dt_est);
	}
      }
    }
  }  // if (GridRank>2)

  // all of the units required to put this into the appropriate enzo internal
  // units, scaled correctly. Note that this does NOT contain a factor
  // of 1/mu, since we don't necessarily know anything about the gas in question.
  all_units = POW(dx,2.0)*POW(LengthUnits,2.0)*DensityUnits*kboltz /
    ( 6.0e-7 * SpitzerFraction * mh * TimeUnits );
  
  dt *= float(all_units);
  dt *= ConductionCourantSafetyNumber;  // for stability, this has to be < 0.5

  if (SpeedOfLightTimeStepLimit) {
    light_cross_time = dx * VelocityUnits / clight;
    dt = max(dt, light_cross_time);
  }

  delete [] Temp;

  return SUCCESS;
 
}
