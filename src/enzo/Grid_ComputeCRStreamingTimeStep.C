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
#include "hydro_rk/EOS.h"



// Member functions
int grid::ComputeCRStreamingTimeStep (float &dt) {
  if (ProcessorNumber != MyProcessorNumber) {return SUCCESS;}
  if (NumberOfBaryonFields == 0) {return SUCCESS;}
  this->DebugCheck("ComputeCRDiffusionTimeStep");

  // Some locals
  float eint, p, rho, v_stream, v_gas, h, dpdrho, dpde, B2,dt_est, hc;

  int size = 1, idx; 
  for (int dim = 0; dim < GridRank; dim++) {size *= GridDimension[dim];};

  FLOAT dx = CellWidth[0][0];

  // We obtain the current cr field ...                                                                                                
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, CRNum, B1Num, B2Num, B3Num, PhiNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum,
				       B1Num, B2Num,B3Num, PhiNum, CRNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Find shortest time scale on the grid patch
  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 1;
    GridEnd[dim] = GridDimension[dim]-2;}

  dt = huge_number;

  for (int k = GridStart[2]; k <= GridEnd[2]; k++)
    for (int j = GridStart[1]; j <= GridEnd[1]; j++)
      for (int i = GridStart[0]; i <= GridEnd[0]; i++){
       
	idx = ELT(i,j,k);
	rho = BaryonField[DensNum][idx];
	
	B2 = BaryonField[B1Num][idx]*BaryonField[B1Num][idx] 
	    + BaryonField[B2Num][idx]*BaryonField[B2Num][idx]
	    + BaryonField[B3Num][idx]*BaryonField[B3Num][idx];
	v_stream = CRStreamVelocityFactor*sqrt(B2/rho);
	hc = CRStreamStabilityFactor * dx; // hard coded
	dt_est = 0.25 * dx*dx / v_stream / hc;
	dt = min(dt, dt_est);

      } // end triple for loop
  return SUCCESS; 
}
