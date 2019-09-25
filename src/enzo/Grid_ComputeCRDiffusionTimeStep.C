/***********************************************************************
/
/  GRID CLASS (Compute Cosmic Ray Diffusion Time-scale)
/
/  written by:  Munier A. Salem
/  date:        January, 2011
/
/  PURPOSE:  Calculates the shortest time scale for cosmic ray diffusion
/            on the current grid.
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
int grid::ComputeCRDiffusionTimeStep (float &dt) {
  if (ProcessorNumber != MyProcessorNumber) {return SUCCESS;}
  if (NumberOfBaryonFields == 0) {return SUCCESS;}
  this->DebugCheck("ComputeCRDiffusionTimeStep");

  // Some locals
  float kappa, dt_est;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  double units = ((double)LengthUnits)*LengthUnits/TimeUnits;

  int size = 1; 
  for (int dim = 0; dim < GridRank; dim++) {size *= GridDimension[dim];};

  FLOAT dx = CellWidth[0][0];


  // Find shortest time scale on the grid patch
  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 1;
    GridEnd[dim] = GridDimension[dim]-2;}

  dt = huge_number;


  /* timestep is calculated as dt < 0.5 * dx^2 / kappa, where
     kappa is the cosmic ray diffusion constant.  */ 
 
  for (int k = GridStart[2]; k <= GridEnd[2]; k++)
    for (int j = GridStart[1]; j <= GridEnd[1]; j++)
      for (int i = GridStart[0]; i <= GridEnd[0]; i++){
	  
  	if( 1 == CRDiffusion ) 	// Constant kappa model
	    kappa = CRkappa/units;

	dt_est = .5 * dx*dx / kappa;

	dt = min(dt, dt_est);

      } // end triple for loop
  return SUCCESS; 
}
