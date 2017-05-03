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
int QuantumGetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

// Member functions
int grid::ComputeQuantumTimeStep (float &dt) {
  if (ProcessorNumber != MyProcessorNumber) {return SUCCESS;}
  if (NumberOfBaryonFields == 0) {return SUCCESS;}
  this->DebugCheck("ComputeQuantumTimeStep");

  // Some locals
  int DensNum, TENum, GENum, Vel1Num, Vel2Num, Vel3Num;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0, aUnits = 1.0;
  FLOAT a = 1.0, dadt;
  double MassUnits = 1.0;
  float hmcoef=1.0;

  // Get system of units
  if (QuantumGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
   hmcoef = 5.9157166856e27*TimeUnits/pow(LengthUnits,2)/FDMMass;


  dt = huge_number;

  /* timestep is calculated as dt < 0.5 * (a*dx)^2 / alpha, where
     alpha = thermal diffusivity = Kappa/(number density * k_boltz).
     The 0.5 is actually put in ComputeTimeStep 
     (as ConductionCourantSafetyFactor), and the rest is calculated
     here.  We first calculate the density and temperature parts in
     the loop (this is the only stuff that varies), and then multiply
     by all of the constants necessary to put this into Enzo internal
     units.  Since the saturation coefficient requires the temperature 
     gradients, this calculation now looks a lot like what is done in 
     Grid_ComputeHeat. */
  
  if (ComovingCoordinates){
      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    ENZO_FAIL("Error in ComputeExpansionFactor.\n");
     }
  }


  FLOAT dx = CellWidth[0][0];

  if (GridRank>1) {
    dx = min( dx, CellWidth[1][0]);
  } // if (GridRank>1)
  
  if (GridRank>2) {
    dx = min( dx, CellWidth[2][0]);
	  } // if (GridRank>2)

  
  dt = pow(dx*a,2)/(4.*hmcoef);
  fprintf(stderr, "time coefficients %lf %f %f %f\n", hmcoef,dx,dt,a);


  return SUCCESS;
 
}
