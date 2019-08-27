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
  int size = 1, idx, i,j,k, Nsub=0; 
  float *cr, *Bx, *By, *Bz, crOld, hc, pcr, rho, dpdrho, dpde, p, h, cs, eint, B2;
  float dtSubcycle, dtSoFar, dCRdt, dCRdt_alt, dCRdx, dCRdy, dCRdz, oldCRdt;
  float bx, by, bz, va_x, va_y, va_z, vs_x, vs_y, vs_z, BdotCR; 

  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *Fcx  = new float[size];
  float *Fcy = new float[size];
  float *Fcz = new float[size];

  float *Fcx_alt  = new float[size];
  float *Fcy_alt = new float[size];
  float *Fcz_alt = new float[size];

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

  // Sub-cycle, computing and applying diffusion

  // dtSubcycle = timestep of this subcycle
  // dtSoFar = overall timestep taken
  // dtFixed = fixed timestep for the entire level.

  dtSoFar = 0.0;

  Nsub=0;  // number of subcycles

  while(dtSoFar < dtFixed){

    // compute this subcycle timestep

    if (this->ComputeCRStreamingTimeStep(dtSubcycle) == FAIL) {
      ENZO_FAIL("Error in ComputeCRDiffusionTimeStep.");
    }
    dtSubcycle *= CRCourantSafetyNumber;  // for stability

    // make sure we don't extend past dtFixed

    dtSubcycle = dtFixed;//min(dtSubcycle, dtFixed-dtSoFar);

    // compute dCR/dt for each cell.

    int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

    /* Set up start and end indexes to cover all of grid except outermost cells. */

    for (int dim = 0; dim<GridRank; dim++ ) {
      GridStart[dim] = 1;
      GridEnd[dim] = GridDimension[dim]-1;
    }

    hc = 0.007; // value taken from Ruzskowski+ 2017

    /* Compute CR fluxes at each cell face. */
      for (k = GridStart[2]; k <= GridEnd[2]; k++)
	for (j = GridStart[1]; j <= GridEnd[1]; j++)
	  for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	    idx = ELT(i,j,k);

	    pcr = (CRgamma - 1.0) *cr[idx];
	    rho = BaryonField[DensNum][idx];

	    
	    dCRdx = (cr[idx]-cr[ELT(i-1,j,k)])/dx[0];
	    dCRdy = (cr[idx]-cr[ELT(i,j-1,k)])/dx[1];
	    dCRdz = (cr[idx]-cr[ELT(i,j,k-1)])/dx[2];

	    // If anisotropic CR streaming, CRs stream with a velocity proportional                 
	    // to the Alfven velocity     
	     
	    va_x = Bx[idx] / sqrt(rho);
	    va_y = By[idx] / sqrt(rho);
	    va_z = Bz[idx] / sqrt(rho);

	    B2 = Bx[idx]*Bx[idx] + By[idx]*By[idx] + Bz[idx]*Bz[idx]; 
	    bx = Bx[idx] / sqrt(B2);
	    by = By[idx] / sqrt(B2);
	    bz = Bz[idx] / sqrt(B2);
	      
	    // Need this to find the sign of B dot dCR/dx
            BdotCR = bx*dCRdx;
	    if (GridRank > 1) BdotCR += by*dCRdy;
	    if (GridRank > 1) BdotCR += bz*dCRdz;
 
	    vs_x = -sign(BdotCR)*va_x;  
	    Fcx[idx] = CRgamma*pcr*vs_x;
	    Fcx_alt[idx] = -CRgamma*pcr*va_x*tanh(0.007*BdotCR / ((cr[idx] == 0.0) ? 1.0 : cr[idx])); 
	    
	    if( GridRank > 1){
	      vs_y = -sign(BdotCR)*va_y; 
	      Fcy[idx] = CRgamma*pcr*vs_y; 
	      Fcy_alt[idx] = -CRgamma*pcr*va_y*tanh(0.007*BdotCR / ((cr[idx] == 0.0) ? 1.0 : cr[idx]));
	    }
	    if( GridRank > 2){
	      vs_z = -sign(BdotCR)*va_z; 
	      Fcz[idx] = CRgamma*pcr*vs_z; 
	      Fcz_alt[idx] = -CRgamma*pcr*va_z*tanh(0.007*BdotCR / ((cr[idx] == 0.0) ? 1.0 : cr[idx]));
	    }

	  } // end triple for

      /* Trim GridEnd so that we don't apply fluxes to cells that don't have
	 them computed on both faces. */

       for (int dim = 0; dim<GridRank; dim++) {
	GridEnd[dim]--;
	}

      // And then update the current CR baryon field (Could be combined with step above)

      for (k = GridStart[2]; k <= GridEnd[2]; k++) 
	for (j = GridStart[1]; j <= GridEnd[1]; j++) 
	  for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	    idx = ELT(i,j,k);
	    crOld = cr[idx];

	    //	    dCRdt = dtSubcycle*(Fcx[ELT(i+1,j,k)]-Fcx[idx])/dx[0];
	    dCRdt = 0.5*dtSubcycle*(Fcx[ELT(i+1,j,k)]-Fcx[ELT(i-1,j,k)])/dx[0]; 
            dCRdt_alt = 0.5*dtSubcycle*(Fcx_alt[ELT(i+1,j,k)]-Fcx_alt[ELT(i-1,j,k)])/dx[0];

	    if( GridRank > 1 ){
	      dCRdt += 0.5*dtSubcycle*(Fcy[ELT(i,j+1,k)]-Fcy[ELT(i,j-1,k)])/dx[1];
	      dCRdt_alt += 0.5*dtSubcycle*(Fcy_alt[ELT(i,j+1,k)]-Fcy_alt[ELT(i,j-1,k)])/dx[1];
	    }
	    if( GridRank > 2 ){
	      dCRdt += 0.5*dtSubcycle*(Fcz[ELT(i,j,k+1)]-Fcz[ELT(i,j,k-1)])/dx[2];
	      dCRdt_alt += 0.5*dtSubcycle*(Fcz_alt[ELT(i,j,k+1)]-Fcz_alt[ELT(i,j,k-1)])/dx[2];
	    }

	    if(((cr[idx] - dCRdt)< 0) || isnan(dCRdt)){
	      oldCRdt = dCRdt;
	      dCRdt = dCRdt_alt;
	      //	      printf("using dCRdt_alt = %"ESYM", oldCRdt = %"ESYM"\n", dCRdt_alt, dCRdt);
     	    } // end err if

	    BaryonField[CRNum][idx] -= CRStreamingFactor*dCRdt;
	    if((cr[idx]< 0) || isnan(cr[idx])){
		     cr[idx] = 0.0; 
	    }

	    } // triple for loop
    // increment timestep
    dtSoFar += dtSubcycle;
    Nsub++;

  } // while(dtSoFar < dtFixed)

  if (debug) 
    printf("Grid::ComputeCRStreaming:  Nsubcycles = %"ISYM", dx=%"ESYM"\n", Nsub, dx[0]); 
	
  delete [] Fcx;
  delete [] Fcy;
  delete [] Fcz;

  delete [] Fcx_alt;
  delete [] Fcy_alt;
  delete [] Fcz_alt;

  return SUCCESS;  
}

