/***********************************************************************
/
/  GRID CLASS (Compute de/dt due to thermal conduction)
/
/  written by:  David A. Ventimiglia and Brian O'Shea
/  date:        December, 2009
/
/  PURPOSE:  Calculates the heat flowing into (or out of) the cells of
/  a grid patch due to thermal conduction.  As of early 2011, it also includes
/  anisotropic conduction (though this only works if MHD is also turned on)
/  Note that we have implemented the anisotropic conduction from Parrish & 
/  Stone 2005 -- NOT the one that's actually in the paper, but the one that
/  they have implemented into Athena.  
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

// Function prototypes
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

/* Limiter functions */
static float minmod(float var1, float var2);  /* minmod limiter */
static float mcd(float var1, float var2);     /* monotonized central diff. limiter */

// Member functions
int grid::ComputeHeat (float dedt[]) {

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("ComputeHeat");

  // Some locals
  int DensNum, TENum, GENum, Vel1Num, Vel2Num, Vel3Num;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0, aUnits = 1.0;
  FLOAT a = 1.0, dadt;
  double MassUnits = 1.0;
  float *rho,*Bx,*By,*Bz;
  double kappa_star = 6.0e-7;
  float Bx_face, By_face, Bz_face, Bxhat, Byhat, Bzhat, Bmag;
  float dTx, dTy, dTz;
  int x_offset=0, y_offset=0, z_offset=0;
  float Bhat; // gotta get rid of this later

  int size = 1, grid_index, right_side_index;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  float *Temp = new float[size];
  FLOAT dx = CellWidth[0][0];

  if (AnisotropicConduction){
    // find magnetic fields: no dimensionality check because they all seem to
    // be initialized, always.
    iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
    // make masks (for convenience)
    Bx = BaryonField[iBx];
    By = BaryonField[iBy];
    Bz = BaryonField[iBz];

    /* offsets are because of the size of the finite-difference stencil for calculating
       cross derivatives in heat flux: you need i,j,k +-1, since it's a 3x3(x3) stencil.
       This means that we're ignoring the outermost cell.  Offset is zero for isotropic
       conduction, because it's only local (well, just the i,j,k and +1 cell on each 
       axis. */
    x_offset=1;
    if(GridRank>1) y_offset=1;
    if(GridRank>2) z_offset=1;

  } // if(AnisotropicConduction)

  // Zero-out the de/dt array
  for (int i=0; i<size; i++) 
    dedt[i] = 0.0;
  
  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, 
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }

  // conversion from CGS to Enzo internal units for de/dt
  double units = a * POW(TimeUnits, 3.0) * POW(aUnits, 2.0) / 
    POW(LengthUnits, 4.0) / DensityUnits;

  // for conduction saturation
  double saturation_factor = 4.874e-20 / (DensityUnits * LengthUnits * dx);
                                        // 4.2 * lambda_e * mH
                                        // lambda_e from Jubelgas ea 2004
                                        // mH for converting rho into n_e
                                        // dx for dT/dx

  // Get field identifiers
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, 
				       Vel2Num, Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  // Get temperature, density, and energy fields
  if (this->ComputeTemperatureField(Temp) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  // mask for baryon density
  rho = BaryonField[DensNum];

  // Set up a struct to hold properties defined on cell faces
  struct cellface {float T, dT, kappa, dedt, rho;} l, r, cfzero;

  // zero struct
  cfzero.T = cfzero.dT = cfzero.kappa = cfzero.dedt = cfzero.rho = 0.0;

  /* When computing heat, loop over the whole grid, EXCEPT for the
     first and last cells.*/
  int GridStart[] = {0, 0, 0}, 
    GridEnd[] = {0, 0, 0};

  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  /* for each grid rank (x,y,z) loop over 1D 'pencils' and calculate
     heat transport into/out of each cell.  We assume that the left 
     face is zero for the first cell, and then set left face of cell
     i+1 to right face of cell i, to save some computation.  We have
     to multiply through by various constants, but this is done 
     at the end.

     The x,y,z offsets in the loops are 0 for isotropic, 1 for 
     anisotropic conduction. 

     Note that the GridRank==1 loop is heavily commented, but the 
     GridRank=2,3 (y and z) are not commented very much.  Read the x
     loop to understand what's going on. */
  if (GridRank>0) {
    for (int k = GridStart[2]+z_offset; k <= GridEnd[2]-z_offset; k++) {
      for (int j = GridStart[1]+y_offset; j <= GridEnd[1]-y_offset; j++) {
	r = cfzero;
	for (int i = GridStart[0]+x_offset; i <= GridEnd[0]-x_offset; i++) {
	  l = r;  // left face of this cell is the right face of the previous cell.
	  r = cfzero;  // zero out right face: we're going to set this later.

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i+1,j,k);

	  if( i != GridEnd[0]-x_offset ){

	    // get temperature, temperature gradient on + face of cell
	    // (the 'l' struct has it on the right face)
	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[right_side_index] - Temp[grid_index];

	    // kappa is the spitzer conductivity, which scales as 
	    // the temperature to the 2.5 power
	    r.kappa = kappa_star*POW(r.T, 2.5);
	    // conduction saturation
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    // modify flux based on magnetic field orientation, if we are using
	    // anisotropic conduction
	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);  // magnetic fields including sign and magnitude
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face + tiny_number*tiny_number), 0.5); // just magnitude
	      Bxhat = Bx_face/Bmag;  // unit vectors
	      Byhat = By_face/Bmag;
	      Bzhat = Bz_face/Bmag;

	      /* calculating derivatives for heat flux:

		 q_vector = -X_c bhat * (bhat dot grad T), so

		 q_x = -Kappa * (bx*bx*dT/dx + bx*by*dT/dy + bx*bz*dT/dz)

		 with the derivatives being partials.  It's analogous for the other dimensions.  For x, the 
		 calculation of dT/dx is straightforward, but you need to calculate the dT/dy (and possibly 
		 dT/dz) derivatives on the +x face of the cell.  If you want to do this in a symmetric way
		 it requires using the j+-1 and k+-1 cells.  We use a monotonize central difference limiter
		 for the off-axis temperature derivatives to avoid spurious oscillations (though I've left the
		 simpler code in here as well, for the sake of friendliness toward future generations).
	       */

	      dTx =Temp[ELT(i+1,j,k)]- Temp[ELT(i,j,k)];

	      //dTy = ((Temp[ELT(i,j+1,k)] - Temp[ELT(i,j-1,k)]) + (Temp[ELT(i+1,j+1,k)] - Temp[ELT(i+1,j-1,k)]))/4.0;
	      dTy = mcd( mcd( Temp[ELT(i+1,j+1,k)] - Temp[ELT(i+1,j,k)], Temp[ELT(i+1,j,k)] - Temp[ELT(i+1,j-1,k)] ),
			 mcd( Temp[ELT(i,  j+1,k)] - Temp[ELT(i,  j,k)], Temp[ELT(i,  j,k)] - Temp[ELT(i,  j-1,k)] ) );

	      // calculate energy flux across this face. Note that the negative sign is missing from the expression
	      // for heat flux: this is accounted for later, in Grid::ConductHeat (where deltat*de/dt is added, not
	      // subtracted.  Also, factors of dx and units are applied at the end of the routine.
	      r.dedt +=  r.kappa*AnisotropicConductionSpitzerFraction*(Bxhat*Bxhat*dTx + Bxhat*Byhat*dTy ); 

	      // If anisotropic conduction is turned on, we _know_ that it's at least 2D, but there's no guarantee that
	      // it's 3D.  So, check, and then if we're using 3D calculate the z cross-derivative
	      if(GridRank > 2){

		//dTz =  ((Temp[ELT(i,j,k+1)] - Temp[ELT(i,j,k-1)]) + (Temp[ELT(i+1,j,k+1)] - Temp[ELT(i+1,j,k-1)]))/4.0;
		dTz = mcd( mcd( Temp[ELT(i+1,j,k+1)] - Temp[ELT(i+1,j,k)], Temp[ELT(i+1,j,k)] - Temp[ELT(i+1,j,k-1)] ), 
			   mcd( Temp[ELT(i,  j,k+1)] - Temp[ELT(i,  j,k)], Temp[ELT(i,  j,k)] - Temp[ELT(i,  j,k-1)] ) );

		r.dedt += r.kappa*AnisotropicConductionSpitzerFraction*Bxhat*Bzhat*dTz;
	      }

	    } // if(AnisotropicConduction) 

	    // add in energy transport because of isotropic conduction as well.
	    if(IsotropicConduction){ 
	      // factors of dx and units done later (also missing negative sign accounted for)
	      r.dedt += r.kappa*IsotropicConductionSpitzerFraction*r.dT;  
	    }

	  }

	  // actually calculate de/dt: note that there's nominally a missing negative sign, which accounts for
	  // the missing sign above
	  dedt[grid_index] += (r.dedt - l.dedt)/rho[grid_index];
	}
      }
    }
  } // if (GridRank>0)


  if (GridRank>1) {
    for (int i = GridStart[0]+x_offset; i <= GridEnd[0]-x_offset; i++) {
      for (int k = GridStart[2]+z_offset; k <= GridEnd[2]-z_offset; k++) {
	r = cfzero;
	for (int j = GridStart[1]+y_offset; j <= GridEnd[1]-y_offset; j++) {
	  l = r;
	  r = cfzero;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j+1,k);

	  if( j != GridEnd[1]-y_offset ){

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[right_side_index] - Temp[grid_index];
	    
	    r.kappa = kappa_star*POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face + tiny_number*tiny_number), 0.5);

	      Bxhat = Bx_face/Bmag;
	      Byhat = By_face/Bmag;
	      Bzhat = Bz_face/Bmag;

	      //dTx = ((Temp[ELT(i+1,j,k)] - Temp[ELT(i-1,j,k)]) + (Temp[ELT(i+1,j+1,k)] - Temp[ELT(i-1,j+1,k)]))/4.0;
	      dTx = mcd( mcd( Temp[ELT(i+1,j+1,k)] - Temp[ELT(i,j+1,k)], Temp[ELT(i,j+1,k)] - Temp[ELT(i-1,j+1,k)] ),
			 mcd( Temp[ELT(i+1,j,  k)] - Temp[ELT(i,j,  k)], Temp[ELT(i,j,  k)] - Temp[ELT(i-1,j,  k)] ) );

	      dTy = Temp[ELT(i,j+1,k)]-Temp[ELT(i,j,k)];

	      r.dedt += r.kappa*AnisotropicConductionSpitzerFraction*(Bxhat*Byhat*dTx + Byhat*Byhat*dTy );  

	      if(GridRank > 2){
		//dTz =  ((Temp[ELT(i,j,k+1)] - Temp[ELT(i,j,k-1)]) + (Temp[ELT(i,j+1,k+1)] - Temp[ELT(i,j+1,k-1)]))/4.0;
		dTz = mcd( mcd( Temp[ELT(i,j+1,k+1)] - Temp[ELT(i,j+1,k)], Temp[ELT(i,j+1,k)] - Temp[ELT(i,j+1,k-1)] ), 
			   mcd( Temp[ELT(i,j,  k+1)] - Temp[ELT(i,j,  k)], Temp[ELT(i,j,  k)] - Temp[ELT(i,j,  k-1)] ) );

		r.dedt += r.kappa*AnisotropicConductionSpitzerFraction*Byhat*Bzhat*dTz;
	      }

	    } // if(AnisotropicConduction)

	    if(IsotropicConduction){
	      r.dedt += r.kappa*IsotropicConductionSpitzerFraction*r.dT;
	    }

	  }

	  dedt[grid_index] += (r.dedt - l.dedt)/rho[grid_index];
	}
      }
    }
  } // if (GridRank>1)
  

  if (GridRank>2) {
    for (int j = GridStart[1]+y_offset; j <= GridEnd[1]-y_offset; j++) {
      for (int i = GridStart[0]+x_offset; i <= GridEnd[0]-x_offset; i++) {
	r = cfzero;
	for (int k = GridStart[2]+z_offset; k <= GridEnd[2]-z_offset; k++) {
	  l = r;
	  r = cfzero;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j,k+1);

	  if( k != GridEnd[2]-z_offset ){

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[right_side_index] - Temp[grid_index];

	    r.kappa = kappa_star*POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face + tiny_number*tiny_number), 0.5);

	      Bxhat = Bx_face/Bmag;
	      Byhat = By_face/Bmag;
	      Bzhat = Bz_face/Bmag;

	      // don't need to check if GridRank>2, because that's the only reason this loop gets called.
	      //dTx = ((Temp[ELT(i+1,j,k)] - Temp[ELT(i-1,j,k)]) + (Temp[ELT(i+1,j,k+1)] - Temp[ELT(i-1,j,k+1)]))/4.0;
	      dTx = mcd( mcd( Temp[ELT(i+1,j,k+1)] - Temp[ELT(i,j,k+1)], Temp[ELT(i,j,k+1)] - Temp[ELT(i-1,j,k+1)] ), 
			 mcd( Temp[ELT(i+1,j,k  )] - Temp[ELT(i,j,k  )], Temp[ELT(i,j,k  )] - Temp[ELT(i-1,j,k  )] ) );

	      //dTy = ((Temp[ELT(i,j+1,k)] - Temp[ELT(i,j-1,k)]) + (Temp[ELT(i,j+1,k+1)] - Temp[ELT(i,j-1,k+1)]))/4.0;
	      dTy = mcd( mcd( Temp[ELT(i,j+1,k+1)] - Temp[ELT(i,j,k+1)], Temp[ELT(i,j,k+1)] - Temp[ELT(i,j-1,k+1)] ), 
			 mcd( Temp[ELT(i,j+1,k  )] - Temp[ELT(i,j,k  )], Temp[ELT(i,j,k  )] - Temp[ELT(i,j-1,k  )] ) );

	      dTz = Temp[ELT(i,j,k+1)]-Temp[ELT(i,j,k)];

	      r.dedt += r.kappa*AnisotropicConductionSpitzerFraction*(Bzhat*Bxhat*dTx + Bzhat*Byhat*dTy + Bzhat*Bzhat*dTz); 

	    } // if(AnisotropicConduction)

	    if(IsotropicConduction){
	      r.dedt += r.kappa*IsotropicConductionSpitzerFraction*r.dT;
	    }

	  }

	  dedt[grid_index] += (r.dedt - l.dedt)/rho[grid_index];
	}
      }
    }
  }  // if (GridRank>2)

  // sweep through once and multiply everything by constants 
  // (units and factor of 1/dx^2)
  for (int i=0; i<size; i++) 
    dedt[i] *= units/(dx*dx);

  // cleanup and return
  delete [] Temp;

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
