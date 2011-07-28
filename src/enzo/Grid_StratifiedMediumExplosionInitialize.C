////////////////////////////////////////////////////////////////////////////////
//
//  GRID CLASS
//
//  written by: Brian O'Shea
//  date:       March 2010
//  modified1:  
//
//  PURPOSE: Initialize conduction bubble test problem.  
//
//  RETURNS: FAIL or SUCCESS
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// Function prototypes
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int FindField(int field, int farray[], int numfields);


// Grid Initializer: all input values are in Enzo internal units
int grid::StratifiedMediumExplosionInitialize(FLOAT BubbleRadius, int PulseType,
				      float ExplosionEnergy, FLOAT BubbleCenter[MAX_DIMENSION]) {

  if (debug) {
    printf("Entering ConductionBubbleInitialize\n");
    fflush(stdout);
  }

  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;

  FLOAT x,y,z, r2;

  int i,j,k;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;

  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, 
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  float delta, this_delta, g;
  
  FLOAT ScaleHeight = 8.0e+5;  // scale height in cm
  float GroundTemp = 300.0;   // ground temperature in K
  float GroundEnergy, GroundDensity;

  float Boltzmann = 1.38e-16, mu = 1.2, mh=1.67e-24;

  GroundEnergy = (Boltzmann*GroundTemp)/((Gamma - 1.0)*mu*mh);
  GroundEnergy /= (VelocityUnits*VelocityUnits);


  GroundDensity = 1.2e-3 / DensityUnits;  // grams/cc in Enzo internal units
  ScaleHeight /= LengthUnits;  // scale height in enzo internal units

  delta = 1.0;  //POW(DeltaEntropy, 0.6);


  // ExplosionEnergy comes in as kilotons
  // one kiloton = 4.184e12 Joules = 4.184e19 ergs
  ExplosionEnergy *= 4.184e19;  // now in ergs

  printf("ExplosionEnergy is %e (pre-conversion)\n",ExplosionEnergy);
  
  // now in Enzo internal units
  ExplosionEnergy /= (double(DensityUnits) * POW(double(LengthUnits),5.0) * POW(double(TimeUnits),-2.0) );
 
  printf("ExplosionEnergy is %e (post-conversion)\n",ExplosionEnergy);

  float Mass;

  if(GridRank==2){
    Mass = GroundDensity * 3.1415 * POW(BubbleRadius,3.0);
  }

  if(GridRank==3){
    Mass = GroundDensity * 4.0*3.1415/3.0 * POW(BubbleRadius,3.0);
  }

  // energy per mass
  ExplosionEnergy /= Mass;

  printf("ExplosionEnergy is %e (mass is %e, bubble radius is %e)\n",ExplosionEnergy, Mass, BubbleRadius);

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  printf("GroundEnergy is %e in code units (explosion energy is %e)\n",GroundEnergy,ExplosionEnergy);

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;


  // calculate gravitational constant in Enzo internal units.
  g = fabs(UniformGravityConstant)*LengthUnits/(TimeUnits*TimeUnits);

  if(UniformGravity==0) g = 0.0;  // if gravity is off make sure it's zero

  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  int ii, small_index;
  FLOAT smallest_d, celldist;

  // loop over grid and set cell values: we're setting both the 
  // pulse values and the background values here.
  for (k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {

	/* Compute position of the cell and find distance 
	   from the center of the bubble */
	x=y=z=0.0;

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	r2 = POW(x-BubbleCenter[0], 2.0);  // currently distance^2

	if(GridRank>1){
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  r2 += POW(y-BubbleCenter[1], 2.0);
	}

	if(GridRank>2){
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  r2 += POW(z-BubbleCenter[2], 2.0);
	}

	celldist = POW(r2,0.5);  // go from distance^2 to just plain distance

	// set baryon density
	BaryonField[DensNum][ELT(i,j,k)] = GroundDensity * exp(-1.0*y/ScaleHeight);  // background

	if(HydroMethod==Zeus_Hydro){  // ZEUS
	  BaryonField[TENum][ELT(i,j,k)] = GroundEnergy ;  // TE = gas energy

	  if(celldist < BubbleRadius) 
	    BaryonField[TENum][ELT(i,j,k)] += ExplosionEnergy;

	} else{ // PPM
	  
	  BaryonField[TENum][ELT(i,j,k)] = GroundEnergy ;

	  if(celldist < BubbleRadius) 
	    BaryonField[TENum][ELT(i,j,k)] += ExplosionEnergy;

	  if(DualEnergyFormalism){

	    BaryonField[GENum][ELT(i,j,k)] = GroundEnergy ;

	    if(celldist < BubbleRadius)     
	      BaryonField[GENum][ELT(i,j,k)] += ExplosionEnergy;

	  } // DEF
	} // PPM

	/* HERE WOULD BE AN APPROPRIATE PLACE TO MODIFY THE VELOCITY FIELD IF YOU NEED TO */

	/* we've messed with the total baryon energy above, so the values are no longer consistent with what was set in
	   Grid::InitializeUniformGrid.  Regardless of whether we've fiddled with the velocity or the magnetic fields, 
	   we need to add that energy back in! */

	// total energy needs to be updated to take into account gas velocity if hydro is PPM or MHD
	// in ZEUS 'total energy' is really internal energy, so we don't have to worry about this.
	if(HydroMethod != Zeus_Hydro){

	  BaryonField[TENum][ELT(i,j,k)] += 0.5*POW(BaryonField[Vel1Num][ELT(i,j,k)], 2.0);
	  if(GridRank > 1)
	    BaryonField[TENum][ELT(i,j,k)] += 0.5*POW(BaryonField[Vel2Num][ELT(i,j,k)], 2.0);
	  if(GridRank > 2)
	    BaryonField[TENum][ELT(i,j,k)] += 0.5*POW(BaryonField[Vel3Num][ELT(i,j,k)], 2.0);
	  
	} // if(HydroMethod != Zeus_Hydro)

	// metallicity
	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE){
	  if(celldist <= BubbleRadius){
	    BaryonField[MetalNum][ELT(i,j,k)] = 
	      BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;
	  } else {
	    BaryonField[MetalNum][ELT(i,j,k)] = tiny_number;
	  }

	  if(j%32==0){  // set lines of tracer field at constant altitude
	    BaryonField[MetalNum][ELT(i,j,k)] = 
	      BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;
	  }
	} // if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)

      } // for(i...)  (loop over all cells in this grid)

  if (debug) {
    printf("Exiting ConductionBubbleInitialize\n");
  }

  return SUCCESS;
}
