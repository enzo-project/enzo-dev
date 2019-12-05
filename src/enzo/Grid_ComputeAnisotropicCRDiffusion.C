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
  int size = 1, idx, i,j,k, Nsub=0; 
  float *cr, *Bx, *By, *Bz,B2, crOld, kappa;
  float dtSubcycle, dtSoFar, dCRdt, dCRdt_tan, dCRdt_iso;
  float dCRx, dCRy, dCRz, dbCRx, dbCRy, dbCRz;
  float Bx_xface, By_xface, Bz_xface, B2_xface;
  float Bx_yface, By_yface, Bz_yface, B2_yface;
  float Bx_zface, By_zface, Bz_zface, B2_zface;

  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];
  float *kdCR = new float[size];
  float *kdCRdx  = new float[size];
  float *kdCRdy = new float[size];
  float *kdCRdz = new float[size];
  float *kdCR_tan = new float[size];
  float *bx = new float[size];
  float *by = new float[size];
  float *bz = new float[size];
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

    kappa = CRkappa/units;        // Constant Kappa Model  

    /* Compute CR fluxes at each cell face. */
    //  for (int NC = 0; NC < 2; NC++){
      for (k = GridStart[2]; k <= GridEnd[2]; k++)
	for (j = GridStart[1]; j <= GridEnd[1]; j++)
	  for (i = GridStart[0]; i <= GridEnd[0]; i++) {
	    idx = ELT(i,j,k);

        
	    Bx_xface = 0.5*(Bx[idx] + Bx[ELT(i-1,j,k)]);
            By_xface = 0.5*(By[idx] + By[ELT(i-1,j,k)]);
	    Bz_xface = 0.5*(Bz[idx] + Bz[ELT(i-1,j,k)]);
	    
	    Bx_yface = 0.5*(Bx[idx] + Bx[ELT(i,j-1,k)]);
	    By_yface = 0.5*(By[idx] + By[ELT(i,j-1,k)]);
	    Bz_yface = 0.5*(Bz[idx] + Bz[ELT(i,j-1,k)]);

	    Bx_zface = 0.5*(Bx[idx] + Bx[ELT(i,j,k-1)]);
	    By_zface = 0.5*(By[idx] + By[ELT(i,j,k-1)]);
	    Bz_zface = 0.5*(Bz[idx] + Bz[ELT(i,j,k-1)]);

	    B2_xface = Bx_xface*Bx_xface + By_xface*By_xface + Bz_xface*Bz_xface;
	    B2_yface = Bx_yface*Bx_yface + By_yface*By_yface + Bz_yface*Bz_yface;
	    B2_zface = Bx_zface*Bx_zface + By_zface*By_zface + Bz_zface*Bz_zface;

	    B2 = Bx[idx]*Bx[idx] + By[idx]*By[idx] + Bz[idx]*Bz[idx];
	    bx[idx] = Bx[idx]/sqrt(B2);
	    by[idx] = By[idx]/sqrt(B2);
	    bz[idx] = Bz[idx]/sqrt(B2);

	    // kdCRdx, kdCRdy, dkCRdz will be used in isotropic case if B = 0 below
	    // kdCR is for the anisotropic case
	    kdCRdx[idx] = kappa*(cr[idx]-cr[ELT(i-1,j,k)])/dx[0]; 
	    kdCR[idx] = (Bx_xface/sqrt(B2_xface))*kdCRdx[idx];
	    //	    kdCR_tan[idx] = Bx_xface*Bx_xface*kdCRdx[idx]/B2_xface;

	    dCRx =  mcd(cr[ELT(i+1,j,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i-1,j,  k)]);
	    kdCR_tan[idx] = kappa*(Bx_xface*Bx_xface /B2_xface)*dCRx/dx[0];

	    if( GridRank > 1 ){
	      kdCRdy[idx] = kappa*(cr[idx]-cr[ELT(i,j-1,k)])/dx[1];
	      kdCR[idx] += (By_yface/sqrt(B2_yface))*kdCRdy[idx];
	      //	      kdCR_tan[idx] += By_yface*By_yface*kdCRdy[idx] / B2_yface;
	 
	      dCRy =  mcd(cr[ELT(i,j+1,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i,j-1,  k)]);
	      kdCR_tan[idx] += kappa*(By_yface*By_yface /B2_yface)*dCRy/dx[1];
	      // y gradient on x face
              if (j < GridEnd[1])
		dCRy = mcd( mcd( cr[ELT(i,  j+1,k)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i,  j-1,k)] ),
			    mcd( cr[ELT(i-1,j+1,k)] - cr[ELT(i-1,j,k)], cr[ELT(i-1,j,k)] - cr[ELT(i-1,j-1,k)] ));
	      else  dCRy = mcd(cr[ELT(i,  j,k)] - cr[ELT(i,  j-1,k)], cr[ELT(i-1,j,k)] - cr[ELT(i-1,j-1,k)]);

	      kdCR_tan[idx] += kappa*(Bx_xface*By_xface /B2_xface)*dCRy/dx[0];
	      
	      
	      // x gradient on yface
	      if (i < GridEnd[0])
		dCRx = mcd( mcd( cr[ELT(i+1,j,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i-1,j,  k)] ),
			    mcd( cr[ELT(i+1,j-1,k)] - cr[ELT(i,j-1,k)], cr[ELT(i,j-1,k)] - cr[ELT(i-1,j-1,k)] ));

	      else  dCRx = mcd(cr[ELT(i,j,  k)] - cr[ELT(i-1,j,  k)], cr[ELT(i,j-1,k)] - cr[ELT(i-1,j-1,k)]);
	      kdCR_tan[idx] += kappa*(By_yface*Bx_yface / B2_yface)*dCRx/dx[1];
	      
	      
	    }
	    if( GridRank > 2 ){
	      kdCRdz[idx] = kappa*(cr[idx]-cr[ELT(i,j,k-1)])/dx[2];
	      kdCR[idx] +=  (Bz_zface/sqrt(B2_zface))*kdCRdz[idx];
	      //	      kdCR_tan[idx] += Bz_zface*Bz_zface*kdCRdz[idx] / B2_zface; 
	      dCRz =  mcd(cr[ELT(i,j,k+1)] - cr[ELT(i,j,k)], cr[ELT(i,j,k)] - cr[ELT(i,j,k-1)]);
	      kdCR_tan[idx] += kappa*(Bz_zface*Bz_zface /B2_zface)*dCRx/dx[2];

	      // x gradient on z face
              if (i < GridEnd[0])
		dCRx = mcd( mcd( cr[ELT(i+1,j,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i-1,j,  k)] ),
			    mcd( cr[ELT(i+1,j,k-1)] - cr[ELT(i,j,k-1)], cr[ELT(i,j,k-1)] - cr[ELT(i-1,j,k-1)] ));
	      else  dCRx = mcd(cr[ELT(i,j, k)] - cr[ELT(i-1,j,  k)], cr[ELT(i,j,k-1)] - cr[ELT(i-1,j,k-1)]);
	      kdCR_tan[idx] += kappa*(Bz_zface*Bx_zface / B2_zface)*dCRx/dx[2];
	      

	      // y gradient on z face
              if (j < GridEnd[1])
		dCRy = mcd( mcd( cr[ELT(i,j+1,  k)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i,j-1,  k)] ),
			    mcd( cr[ELT(i,j+1,k-1)] - cr[ELT(i,j,k-1)], cr[ELT(i,j,k-1)] - cr[ELT(i,j-1,k-1)] ));
	      else dCRy = mcd(cr[ELT(i,j,  k)] - cr[ELT(i,j-1,  k)], cr[ELT(i,j,k-1)] - cr[ELT(i,j-1,k-1)]); 
	      kdCR_tan[idx] += kappa*(Bz_zface*By_zface / B2_zface)*dCRy/dx[2];
	    
	      
	      if (k < GridEnd[2]){
		// z gradient on x face
		dCRz = mcd( mcd( cr[ELT(i,  j,k+1)] - cr[ELT(i,  j,k)], cr[ELT(i,  j,k)] - cr[ELT(i,  j,k-1)] ),
			    mcd( cr[ELT(i-1,j,k+1)] - cr[ELT(i-1,j,k)], cr[ELT(i-1,j,k)] - cr[ELT(i-1,j,k-1)] ));
		kdCR_tan[idx] += kappa*(Bx_xface*Bz_xface / B2_xface)*dCRz/dx[0];

		// z gradient on y face
		dCRz = mcd( mcd( cr[ELT(i,j,  k+1)] - cr[ELT(i,j,  k)], cr[ELT(i,j,  k)] - cr[ELT(i,j,  k-1)] ),
			    mcd( cr[ELT(i,j-1,k+1)] - cr[ELT(i,j-1,k)], cr[ELT(i,j-1,k)] - cr[ELT(i,j-1,k-1)] ));
		kdCR_tan[idx] += kappa*(By_yface*Bz_yface / B2_yface)*dCRz/dx[1];
	      }

	      else{
                dCRz = mcd(cr[ELT(i,  j,k)] - cr[ELT(i,  j,k-1)],cr[ELT(i-1,j,k)] - cr[ELT(i-1,j,k-1)]);
                kdCR_tan[idx] += kappa*(Bx_xface*Bz_xface / B2_xface)*dCRz/dx[0];

                dCRz = mcd(cr[ELT(i,j,  k)] - cr[ELT(i,j,  k-1)], cr[ELT(i,j-1,k)] - cr[ELT(i,j-1,k-1)]);
                kdCR_tan[idx] += kappa*(By_yface*Bz_yface / B2_yface)*dCRz/dx[1];
	      }
	    } // end GridRank > 2

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

	    B2 = Bx[idx]*Bx[idx] + By[idx]*By[idx] + Bz[idx]*Bz[idx];
	    // if negligibly weak or zero B-field, approximate isotropic diffusion
	    if (B2 == 0.0){
	      dCRdt = dtSubcycle*(kdCRdx[ELT(i+1,j,k)]-kdCRdx[idx])/dx[0];
	      if( GridRank > 1 )
		dCRdt += dtSubcycle*(kdCRdy[ELT(i,j+1,k)]-kdCRdy[idx])/dx[1];
	      if( GridRank > 2 )
		dCRdt += dtSubcycle*(kdCRdz[ELT(i,j,k+1)]-kdCRdz[idx])/dx[2];
	    }
	    
	    // else, do anisotropic diffusion
	    else{
	      dCRdt = dtSubcycle*(bx[ELT(i+1, j, k)]*kdCR[ELT(i+1, j, k)] - bx[idx]*kdCR[idx])/dx[0];
	      if(GridRank > 1)
		dCRdt += dtSubcycle*(by[ELT(i, j+1, k)]*kdCR[ELT(i, j+1, k)] - by[idx]*kdCR[idx])/dx[1];
	      if(GridRank > 2)
		dCRdt += dtSubcycle*(bz[ELT(i, j, k+1)]*kdCR[ELT(i, j, k+1)] - bz[idx]* kdCR[idx])/dx[2];

	      if ((cr[idx] + dCRdt) < 0){ 
		//printf("cr < 0: dCRdt = %"ESYM"    ->", dCRdt);
		
		dCRdt = dtSubcycle*kdCR_tan[idx]/dx[0];
					//printf("New dCRdt = %"ESYM", kdCR_tan[idx] = %"ESYM"\n", dCRdt, kdCR_tan[idx]);
	      }
	     
	    }

	    BaryonField[CRNum][idx] += dCRdt;

            if((cr[idx] < 0) || isnan(cr[idx])){
              printf("CR < 0 (after diff), dtSubcycle = %"ESYM",i,j,k = (%"ISYM", %"ISYM", %"ISYM"), grid lims = \
(%"ISYM", %"ISYM", %"ISYM"), (%"ISYM", %"ISYM", %"ISYM")\n",
		     dtSubcycle, i, j, k, GridStart[0], 
		     GridStart[1], GridStart[2], GridEnd[0], GridEnd[1], GridEnd[2]);	      
	      printf("\t\t>> Old CR: %"ESYM"\n",crOld);
	     
	    } // end err if
	  } // triple for loop
      //} // end NC loop
    // increment timestep
    dtSoFar += dtSubcycle;
    Nsub++;

  } // while(dtSoFar < dtFixed)

  if (debug) 
    printf("Grid::ComputeAnisotropicCRDiffusion:  Nsubcycles = %"ISYM", kappa = %"ESYM", dx=%"ESYM"\n", Nsub, kappa, dx[0]); 
	
  delete [] kdCR;
  delete [] kdCRdx;
  delete [] kdCRdy;
  delete [] kdCRdz;
  delete [] kdCR_tan;
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

