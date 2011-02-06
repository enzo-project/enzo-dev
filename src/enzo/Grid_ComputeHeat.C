/***********************************************************************
/
/  GRID CLASS (Compute de/dt due to thermal conduction)
/
/  written by:  David A. Ventimiglia and Brian O'Shea
/  date:        December, 2009
/
/  PURPOSE:  Calculates the heat flowing into (or out of) the cells of
/  a grid patch due to thermal conduction.
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
  double kappa_star = 6.0e-7 * ConductionSpitzerFraction;
  float Bx_face, By_face, Bz_face, Bhat, Bmag;

  int size = 1, grid_index, right_side_index;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  float *Temp = new float[size];
  FLOAT dx = CellWidth[0][0];

  if (AnisotropicConduction){
    // find fields
    iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
    // make masks (easier later)
    Bx = BaryonField[iBx];
    By = BaryonField[iBy];
    Bz = BaryonField[iBz];

  }

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
     at the end. */
  if (GridRank>0) {
    for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
      for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
	r = cfzero;
	for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i+1,j,k);

	  if(i == GridEnd[0]){
	    r = cfzero;
	  } else {

	    // get temperature, temperature gradient on + face of cell
	    // (the 'l' struct has it on the right face)
	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];

	    // kappa is the spitzer conductivity, which scales as 
	    // the temperature to the 2.5 power
	    r.kappa = kappa_star*POW(r.T, 2.5);
	    // conduction saturation
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    // modify flux based on magnetic field orientation, if we are using
	    // anisotropic conduction
	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face), 0.5);
	      Bhat = fabs(Bx_face)/Bmag;
	      r.kappa *= Bhat;
	    } 

	    r.dedt = r.kappa*r.dT;  // factors of dx and units done later.

	  }

	  dedt[grid_index] += (l.dedt - r.dedt)/rho[grid_index];
	}
      }
    }
  } // if (GridRank>0)

  if (GridRank>1) {
    for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
      for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
	r = cfzero;
	for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j+1,k);

	  if(j==GridEnd[1]){
	    r = cfzero;
	  } else {

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];
	    
	    r.kappa = kappa_star*POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    // modify flux based on magnetic field orientation, if we are using
	    // anisotropic conduction
	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face), 0.5);
	      Bhat = fabs(By_face)/Bmag;
	      r.kappa *= Bhat;
	    } 


	    r.dedt = r.kappa*r.dT;
	  }

	  dedt[grid_index] += (l.dedt - r.dedt)/rho[grid_index];
	}
      }
    }
  } // if (GridRank>1)
  
  if (GridRank>2) {
    for (int j = GridStart[1]; j <= GridEnd[1]; j++) {
      for (int i = GridStart[0]; i <= GridEnd[0]; i++) {
	r = cfzero;
	for (int k = GridStart[2]; k <= GridEnd[2]; k++) {
	  l = r;

	  grid_index = ELT(i,j,k);
	  right_side_index = ELT(i,j,k+1);

	  if(k==GridEnd[2]){
	    r = cfzero;
	  } else {

	    r.T = 0.5 * (Temp[grid_index] + Temp[right_side_index]);
	    r.rho = 0.5 * (rho[grid_index] + rho[right_side_index]);
	    r.dT = Temp[grid_index] - Temp[right_side_index];

	    r.kappa = kappa_star*POW(r.T, 2.5);
	    r.kappa /= (1 + (saturation_factor * r.T * fabs(r.dT) / r.rho));

	    // modify flux based on magnetic field orientation, if we are using
	    // anisotropic conduction
	    if(AnisotropicConduction){
	      Bx_face = 0.5*(Bx[grid_index]+Bx[right_side_index]);
	      By_face = 0.5*(By[grid_index]+By[right_side_index]);
	      Bz_face = 0.5*(Bz[grid_index]+Bz[right_side_index]);
	      Bmag = POW( (Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face), 0.5);
	      Bhat = fabs(Bz_face)/Bmag;
	      r.kappa *= Bhat;
	    } 


	    r.dedt = r.kappa*r.dT;
	  }

	  dedt[grid_index] += (l.dedt - r.dedt)/rho[grid_index];
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
