/***********************************************************************
/
/  GRID CLASS (Compute and apply anisotropic cosmic ray diffusion)
/
/  written by:  Iryna Butsky 
/  date:        March 2017
/
/  PURPOSE:  Calculates and explicit anisotropic cosmic ray diffusion 
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


/* Limiter functions */
static float minmod(float var1, float var2);  /* minmod limiter */
static float mcd(float var1, float var2);     /* monotonized central diff. limiter */

int grid::ComputeAnisotropicCRDiffusion(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k;
  float *cr, *Bx, *By, *Bz, B2, crOld, dCRdt, dCRdt_tan, kappa;
  float dEcrdy_x, dEcrdz_x, dEcrdx_y, dEcrdz_y, dEcrdx_z, dEcrdy_z;
  float bx_xface, by_xface, bz_xface;
  float bx_yface, by_yface, bz_yface;
  float bx_zface, by_zface, bz_zface;

  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];
  float *dEcrdx         = new float[size];
  float *dEcrdy         = new float[size];
  float *dEcrdz         = new float[size];
  float *BdotDelEcr     = new float[size];
  float *BdotDelEcr_tan = new float[size];
  float *bx             = new float[size];
  float *by             = new float[size];
  float *bz             = new float[size];


  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, CRNum, B1Num, B2Num, B3Num, PhiNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum,
					   B1Num, B2Num,B3Num, PhiNum, CRNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  cr = BaryonField[CRNum];
  Bx = BaryonField[B1Num]; 
  By = BaryonField[B2Num]; 
  Bz = BaryonField[B3Num]; 

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
  kappa = CRkappa/units;        // Constant Kappa Model  

  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
  /* Set up start and end indexes to cover all of grid except outermost cells. */
  for (int dim = 0; dim<GridRank; dim++ ) {
    GridStart[dim] = GridStartIndex[dim] - 1;
    GridEnd[dim] = GridEndIndex[dim] + 1;
  }
  /* Compute CR fluxes at each cell face. */
  for (k = GridStart[2]; k <= GridEnd[2]; k++)
    for (j = GridStart[1]; j <= GridEnd[1]; j++)
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	idx = ELT(i,j,k);

	B2 = Bx[idx]*Bx[idx] + By[idx]*By[idx] + Bz[idx]*Bz[idx];
	bx[idx] = Bx[idx] / sqrt(B2);
	by[idx] = By[idx] / sqrt(B2);
	bz[idx] = Bz[idx] / sqrt(B2);

	bx_xface = 0.5*(bx[idx] + bx[ELT(i-1,j,k)]);
        by_xface = 0.5*(by[idx] + by[ELT(i-1,j,k)]);
	bz_xface = 0.5*(bz[idx] + bz[ELT(i-1,j,k)]);
	    
	bx_yface = 0.5*(bx[idx] + bx[ELT(i,j-1,k)]);
        by_yface = 0.5*(by[idx] + by[ELT(i,j-1,k)]);
	bz_yface = 0.5*(bz[idx] + bz[ELT(i,j-1,k)]);

	bx_zface = 0.5*(bx[idx] + bx[ELT(i,j,k-1)]);
	by_zface = 0.5*(by[idx] + by[ELT(i,j,k-1)]);
	bz_zface = 0.5*(bz[idx] + bz[ELT(i,j,k-1)]);


	// CR flux through x-face 
	dEcrdx[idx] = (cr[ELT(i, j, k)] - cr[ELT(i-1, j, k)]) / dx[0];
	if (GridRank > 1){
	  dEcrdy_x = mcd( mcd( cr[ELT(i,  j+1,k)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i,  j-1,k)] ),
                          mcd( cr[ELT(i-1,j+1,k)] - cr[ELT(i-1,j,k)], cr[ELT(i-1,j,k)] - cr[ELT(i-1,j-1,k)] ));
	  dEcrdy_x /= dx[1];
	}
	if (GridRank > 2){
	  dEcrdz_x = mcd( mcd( cr[ELT(i,  j,k+1)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i,  j,k-1)] ),
                          mcd( cr[ELT(i-1,j,k+1)] - cr[ELT(i-1,j,k)], cr[ELT(i-1,j,k)] - cr[ELT(i-1,j,k-1)] ));
	  dEcrdz_x /= dx[2];
	  
	}
	BdotDelEcr[idx] = bx_xface*dEcrdx[idx];
	BdotDelEcr_tan[idx] = by_xface*dEcrdy_x + bz_xface*dEcrdz_x;


	// CR flux through y-face 
	if (GridRank > 1){
	  dEcrdy[idx] = (cr[ELT(i, j, k)] - cr[ELT(i, j-1, k)]) / dx[1];
	  dEcrdx_y =  mcd( mcd( cr[ELT(i+1,  j,k)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i-1,  j,k)] ),
                           mcd( cr[ELT(i+1,j-1,k)] - cr[ELT(i,j-1,k)], cr[ELT(i,j-1,k)] - cr[ELT(i-1,j-1,k)] ));
	  dEcrdx_y /= dx[0];

	  if (GridRank > 2){
	    dEcrdz_y =  mcd( mcd( cr[ELT(i,  j,k+1)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i,  j,k-1)] ),
                             mcd( cr[ELT(i,j-1,k+1)] - cr[ELT(i,j-1,k)], cr[ELT(i,j-1,k)] - cr[ELT(i,j-1,k-1)] ));
	    dEcrdz_y /= dx[2];
	  }
	  BdotDelEcr[idx] += by_yface*dEcrdy[idx];
	  BdotDelEcr_tan[idx] += (bx_yface*dEcrdx_y + bz_yface*dEcrdz_y);
	}

	// CR flux through z-face
	if (GridRank > 2){
	  dEcrdz[idx] = (cr[ELT(i, j, k)] - cr[ELT(i, j, k-1)]) /dx[2];
	  dEcrdx_z =  mcd( mcd( cr[ELT(i+1,j,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i-1,j,  k)] ),
                           mcd( cr[ELT(i+1,j,k-1)] - cr[ELT(i,j,k-1)], cr[ELT(i,j,k-1)] - cr[ELT(i-1,j,k-1)] ));
	  dEcrdx_z /= dx[0];

	  dEcrdy_z = mcd( mcd( cr[ELT(i,j+1,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i,j-1,  k)] ),
                          mcd( cr[ELT(i,j+1,k-1)] - cr[ELT(i,j,k-1)], cr[ELT(i,j,k-1)] - cr[ELT(i,j-1,k-1)] ));
	  dEcrdy_z /= dx[1];

       	  BdotDelEcr[idx] += bz_zface*dEcrdz[idx];
	  BdotDelEcr_tan[idx] += (bx_zface*dEcrdx_z + by_zface*dEcrdy_z);
       }	

      } // end triple for

  for (int dim = 0; dim<GridRank; dim++) {
    GridEnd[dim]--;
  }

  for (k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	idx = ELT(i,j,k);
	crOld = cr[idx];
	
	// if negligibly weak or zero B-field, approximate isotropic diffusion
	if (B2 < tiny_number){
	  dCRdt =    (dEcrdx[ELT(i+1,j,k)]-dEcrdx[idx])/dx[0];
	  if( GridRank > 1 )
	    dCRdt += (dEcrdy[ELT(i,j+1,k)]-dEcrdy[idx])/dx[1];
	  if( GridRank > 2 )
	    dCRdt += (dEcrdz[ELT(i,j,k+1)]-dEcrdz[idx])/dx[2];
	}
	    
	// else, do anisotropic diffusion
	else{
	  dCRdt    = (bx[ELT(i+1, j, k)] * BdotDelEcr[ELT(i+1, j, k)] - bx[idx] * BdotDelEcr[idx]) / dx[0];
	  if(GridRank > 1)
	    dCRdt += (by[ELT(i, j+1, k)] * BdotDelEcr[ELT(i, j+1, k)] - by[idx] * BdotDelEcr[idx]) / dx[1];
	  if(GridRank > 2)
	    dCRdt += (bz[ELT(i, j, k+1)] * BdotDelEcr[ELT(i, j, k+1)] - bz[idx] * BdotDelEcr[idx]) / dx[2];
	}

	BaryonField[CRNum][idx] += kappa * dCRdt * dtFixed;

	// if traditional anisotropic diffusion approach gives unphysical flux, then apply tangential component
	if (BaryonField[CRNum][idx] < 0){
	  dCRdt_tan    = (bx[ELT(i+1, j, k)] * BdotDelEcr_tan[ELT(i+1, j, k)] - bx[idx] * BdotDelEcr_tan[idx]) / dx[0];
	  if(GridRank > 1)
            dCRdt_tan += (by[ELT(i, j+1, k)] * BdotDelEcr_tan[ELT(i, j+1, k)] - by[idx] * BdotDelEcr_tan[idx]) / dx[1];
          if(GridRank > 2)
            dCRdt_tan += (bz[ELT(i, j, k+1)] * BdotDelEcr_tan[ELT(i, j, k+1)] - bz[idx] * BdotDelEcr_tan[idx]) / dx[2];
	  BaryonField[CRNum][idx] += kappa * dCRdt_tan * dtFixed;
	}

        if((cr[idx] < 0) || isnan(cr[idx])){
              printf("CR = %e < 0 (after diff), i,j,k = (%"ISYM", %"ISYM", %"ISYM"), \
                      grid lims = (%"ISYM", %"ISYM", %"ISYM"), (%"ISYM", %"ISYM", %"ISYM")\n",
		     cr[idx], i, j, k, GridStart[0], GridStart[1], GridStart[2], 
		     GridEnd[0], GridEnd[1], GridEnd[2]);	      
	      printf("\t\t>> Old CR: %"ESYM"\n",crOld);
	      cr[idx] = tiny_number; 
	} // end err if
      } // triple for loop

  delete [] dx; 
  delete [] dEcrdx;
  delete [] dEcrdy;
  delete [] dEcrdz;
  delete [] BdotDelEcr;
  delete [] BdotDelEcr_tan;
  delete [] bx;
  delete [] by; 
  delete [] bz;
  return SUCCESS;  
}

/* monotonized central difference limiter - Van Leer 1977. Uses minmod                                                                   
   limiter for convenience. */
static float mcd(float var1, float var2){
  float limited, temp;
  temp = 2.0*minmod(var1, var2);
  limited = minmod(temp,(var1+var2)/2.0);
  return limited;
}

/* minmod limiter (Roe 1986) */
static float minmod(float var1, float var2){
  /*if the variables are of the same sign, return the minimum,                                                                           
   *if the variables are of opposite sign, return zero*/
  float limited;
  if(var1*var2 > 0.0){
    if(fabs(var1) > fabs(var2))
      limited = var2;
    else
      limited = var1;
  }
  else
    limited = 0.0;
  return limited;
}

