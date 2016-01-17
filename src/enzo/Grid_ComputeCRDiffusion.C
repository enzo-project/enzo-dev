/***********************************************************************
/
/  GRID CLASS (Compute and apply cosmic ray diffusion)
/
/  written by:  Munier A. Salem
/  date:        January, 2011
/
/  PURPOSE:  Calculates and applies cosmic ray diffusion 
/  
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
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

int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, double *MassUnits, FLOAT Time);


int grid::ComputeCRDiffusion(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k, Nsub=0; 
  float *cr, crOld, kappa;
  float dtSubcycle, dtSoFar;

  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *dCRdt = new float[size];
  float *kdCRdx  = new float[size];
  float *kdCRdy = new float[size];
  float *kdCRdz = new float[size];

  // We obtain the current cr field ...
	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, CRNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
            Vel3Num, TENum, CRNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  cr = BaryonField[CRNum];

 // Some locals
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }


  double units = ((double)LengthUnits)*LengthUnits/((double)TimeUnits);

  // Sub-cycle, computing and applying diffusion

  // dtSubcycle = timestep of this subcycle
  // dtSoFar = overall timestep taken
  // dtFixed = fixed timestep for the entire level.

  dtSoFar = 0.0;

  Nsub=0;  // number of subcycles

  while(dtSoFar < dtFixed){

    // compute this subcycle timestep

    if (this->ComputeCRDiffusionTimeStep(dtSubcycle) == FAIL) {
      ENZO_FAIL("Error in ComputeConductionTimeStep.");
    }
    dtSubcycle *= CRCourantSafetyNumber;  // for stability

    // make sure we don't extend past dtFixed

    dtSubcycle = min(dtSubcycle, dtFixed-dtSoFar);

    // compute dCR/dt for each cell.

    int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

    /* Set up start and end indexes to cover all of grid except outermost cells. */

    for (int dim = 0; dim<GridRank; dim++ ) {
      GridStart[dim] = 1;
      GridEnd[dim] = GridDimension[dim]-1;
    }

    /* Compute CR fluxes at each cell face. */

    for (k = GridStart[2]; k <= GridEnd[2]; k++)
      for (j = GridStart[1]; j <= GridEnd[1]; j++)
        for (i = GridStart[0]; i <= GridEnd[0]; i++) {
 	  idx = ELT(i,j,k);

	  if( 1 == CRDiffusion )
	    kappa = CRkappa/units;	// Constant Kappa Model

	  kdCRdx[idx] = kappa*(cr[idx]-cr[ELT(i-1,j,k)])/dx[0];
	  if( GridRank > 1 )
	    kdCRdy[idx] = kappa*(cr[idx]-cr[ELT(i,j-1,k)])/dx[1];
	  if( GridRank > 2 )
	    kdCRdz[idx] = kappa*(cr[idx]-cr[ELT(i,j,k-1)])/dx[2];
        } // end triple for

    /* Trim GridEnd so that we don't apply fluxes to cells that don't have
       them computed on both faces. */

    for (int dim = 0; dim<GridRank; dim++) {
      GridEnd[dim]--;
    }

    /* Loop over all all cells and compute cell updats (flux differences) */

    for (k = GridStart[2]; k <= GridEnd[2]; k++)
      for (j = GridStart[1]; j <= GridEnd[1]; j++)
        for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	  idx = ELT(i,j,k);
		
	  dCRdt[idx] = (kdCRdx[ELT(i+1,j,k)]-kdCRdx[idx])/dx[0];
	  if( GridRank > 1 )
	    dCRdt[idx] += (kdCRdy[ELT(i,j+1,k)]-kdCRdy[idx])/dx[1];
	  if( GridRank > 2 )
	    dCRdt[idx] += (kdCRdz[ELT(i,j,k+1)]-kdCRdz[idx])/dx[2];
        }// end triple for

    // And then update the current CR baryon field (Could be combined with step above)

    for (k = GridStart[2]; k <= GridEnd[2]; k++) 
      for (j = GridStart[1]; j <= GridEnd[1]; j++) 
  	for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	  idx = ELT(i,j,k);
	  crOld = cr[idx];
  	  cr[idx] += dCRdt[idx]*dtSubcycle;
          if( cr[idx] < 0.0 ){
            printf("Negative CR (after diff) i,j,k = %"ISYM", %"ISYM", %"ISYM"\n",i,j,k);
            printf("\t\t>> Old CR: %"ESYM"\n",crOld);
          } // end err if
	} // triple for loop

    // increment timestep
    dtSoFar += dtSubcycle;
    Nsub++;

  } // while(dtSoFar < dtFixed)

  if (debug) 
    printf("Grid::ComputeCRDiffusion:  Nsubcycles = %"ISYM", kappa = %"ESYM", dx=%"ESYM"\n", Nsub, kappa, dx[0]); 
	
  delete [] dCRdt;
  delete [] kdCRdx;	
  delete [] kdCRdy;
  delete [] kdCRdz;
  return SUCCESS;  
}

