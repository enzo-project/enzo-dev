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
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

// Member functions
int grid::ComputeConductionTimeStep (float &dt) {
  if (ProcessorNumber != MyProcessorNumber) {return SUCCESS;}
  if (NumberOfBaryonFields == 0) {return SUCCESS;}
  this->DebugCheck("ComputeConductionTimeStep");

  // Some locals
  int DensNum, TENum, GENum, Vel1Num, Vel2Num, Vel3Num;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;
  float dt_est;
  double all_units;

  int size = 1; 
  for (int dim = 0; dim < GridRank; dim++) {size *= GridDimension[dim];};

  float *Temp = new float[size];

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  // get field identifiers
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  // get temperature
  if (this->ComputeTemperatureField(Temp) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  // Find shortest time scale on the grid patch
  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
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
     units.  */
  for (int k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (int j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (int i = GridStart[0]; i <= GridEnd[0]; i++) {

	int idx = ELT(i,j,k);

	dt_est =  BaryonField[DensNum][idx] / POW(Temp[idx], 2.5);

	dt = min(dt, dt_est);
	
      }

  // all of the units required to put this into the appropriate enzo internal
  // units, scaled correctly. Note that this does NOT contain a factor
  // of 1/mu, since we don't necessarily know anything about the gas in question.
  all_units = POW(CellWidth[0][0],2.0)*POW(LengthUnits,2.0)*DensityUnits*kboltz
    / ( 6.0e-7 * ConductionSpitzerFraction * mh * TimeUnits );
  
  dt *= float(all_units);

  delete [] Temp;

  return SUCCESS;
 
}
