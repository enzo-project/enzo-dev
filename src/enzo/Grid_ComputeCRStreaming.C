/***********************************************************************
/
/  GRID CLASS (Compute and apply cosmic ray streaming)
/
/  written by:  Iryna Butsky 
/  date:        August 2017
/
/  PURPOSE:  Calculates and explicit cosmic ray streaming
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


int grid::ComputeCRStreaming(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k; 
  float *cr, *Bx, *By, *Bz, crOld, rho, stream_factor;
  float dCRdt, dCRdx, dCRdy, dCRdz, va_x, va_y, va_z;
  
  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *Fcx = new float[size];
  float *Fcy = new float[size];
  float *Fcz = new float[size];

  // We obtain the current cr field ...
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, CRNum, B1Num, B2Num, B3Num, PhiNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum,
					   B1Num, B2Num,B3Num, PhiNum, CRNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  cr = BaryonField[CRNum];
  Bx = BaryonField[B1Num]; 
  By = BaryonField[B2Num]; 
  Bz = BaryonField[B3Num]; 


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

	rho = BaryonField[DensNum][idx];	    
	dCRdx = (cr[idx]-cr[ELT(i-1,j,k)])/dx[0];
	dCRdy = (cr[idx]-cr[ELT(i,j-1,k)])/dx[1];
	dCRdz = (cr[idx]-cr[ELT(i,j,k-1)])/dx[2];

	// CRs stream with a velocity proportional to the Alfven velocity     
	va_x = Bx[idx] / sqrt(rho);
	va_y = By[idx] / sqrt(rho);
	va_z = Bz[idx] / sqrt(rho);
	
	// Calculate CR flux
	Fcx[idx] = CRgamma * cr[idx] * fabs(va_x) * tanh(dCRdx / CRStreamStabilityFactor);
	if( GridRank > 1)
	  Fcy[idx] = CRgamma * cr[idx] * fabs(va_y) * tanh(dCRdy / CRStreamStabilityFactor);
	
	if( GridRank > 2)
	  Fcz[idx] = CRgamma * cr[idx] * fabs(va_z) * tanh(dCRdz / CRStreamStabilityFactor);

      } // end triple for

  /* Trim GridEnd so that we don't apply fluxes to cells that don't have
     them computed on both faces. */

  for (int dim = 0; dim<GridRank; dim++) {
    GridEnd[dim]--;
  }

  for (k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	idx = ELT(i,j,k);
	crOld = cr[idx];

	dCRdt = 0.5 * (Fcx[ELT(i+1,j,k)]-Fcx[ELT(i-1,j,k)])/dx[0]; 
	if( GridRank > 1 )
	  dCRdt +=  0.5 * (Fcy[ELT(i,j+1,k)]-Fcy[ELT(i,j-1,k)])/dx[1];
	if( GridRank > 2 )
	  dCRdt += 0.5*(Fcz[ELT(i,j,k+1)]-Fcz[ELT(i,j,k-1)])/dx[2];

	BaryonField[CRNum][idx] += CRStreamVelocityFactor * dCRdt * dtFixed;

	if((cr[idx]< 0) || isnan(cr[idx])){
	      printf("CR = %e < 0 (after stream), i,j,k = (%"ISYM", %"ISYM", %"
		     ISYM"), grid lims = (%"ISYM", %"ISYM", %"ISYM"), (%"ISYM", %"ISYM", %"ISYM")\n",
                     cr[idx], i, j, k, GridStart[0],
                     GridStart[1], GridStart[2], GridEnd[0], GridEnd[1], GridEnd[2]);		  
	      cr[idx] = CRdensFloor; 
	}
      } // triple for loop

	
  delete [] Fcx;
  delete [] Fcy;
  delete [] Fcz;

  return SUCCESS;  
}

