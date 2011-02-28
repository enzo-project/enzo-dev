////////////////////////////////////////////////////////////////////////////////
//
//  GRID CLASS
//
//  written by: David A. Ventimiglia & Brian O'Shea
//  date:       March 2010
//  modified1:  
//
//  PURPOSE: 
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

// Grid Initializer
int grid::ConductionTestInitialize (float PulseHeight, FLOAT PulseWidth, int PulseType, FLOAT PulseCenter[MAX_DIMENSION], int FieldGeometry, float BField) {

  if (debug) {
    printf("Entering ConductionTestInitialize\n");
    fflush(stdout);
  }

  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  if(AnisotropicConduction){
    iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
  }

  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  FLOAT sig2 = PulseWidth*PulseWidth;

  FLOAT x,y,z, r2, celldist;

  int i,j,k;

  // loop over grid and set pulse values
  for (k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {

	/* Compute position */
	x=y=z=0.0;

	/* Find distance from center. */

	// radius squared: assume we always want to be at center of 
	// box
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	r2 = POW(x-PulseCenter[0], 2.0);

	if(GridRank>1){
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  r2 += POW(y-PulseCenter[1], 2.0);
	}

	if(GridRank>2){
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  r2 += POW(z-PulseCenter[2], 2.0);
	}

	celldist = POW(r2,0.5);

	float val=1.0;

	/* now we select between our different pulse options */

	if(PulseType == 1){  // gaussian pulse

	  val = 1.0 + exp(-1.0*r2/sig2/2.0)*(PulseHeight-1.0);

	} else if(PulseType == 2){ // square pulse

	  if(r2 <= sig2)
	    val = PulseHeight;
	  
	} else if(PulseType == 3){  // sinusoidal pulse with values along x-axis

	  val = 1.0 + PulseHeight + PulseHeight * sin(2.0 * 3.14158 * x / PulseWidth);
	  
	} else if(PulseType == 4){  // square pulse with smoothed edges (as suggested by A. Kravtsov)

	  val = 1.0 + (PulseHeight-1.0)*(1.0 - tanh((10.*(celldist/PulseWidth-1.0)))) / 2.0;

	} else {

	  ENZO_FAIL("Grid::ConductionTestInitialize: PulseType is not 1,2 or 3!");
	  
	}

	if(HydroMethod==Zeus_Hydro){  // ZEUS
	  BaryonField[TENum][ELT(i,j,k)] *= val;  // TE = gas energy
	} else{ // PPM
	  
	  BaryonField[TENum][ELT(i,j,k)] *= val;  // TE = total energy energy, but velocity=0 here.

	  if(DualEnergyFormalism)
	    BaryonField[GENum][ELT(i,j,k)] *= val;  // if DEF=1, need to separately set the gas internal energy.
	}

	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
	  BaryonField[MetalNum][ELT(i,j,k)] = 
	    BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;

	if(FieldGeometry==1){

	  if(celldist > tiny_number){
	    BaryonField[iBx][ELT(i,j,k)] = -1.0*BField*(y-.5)/POW( (POW(x-.5,2.0)+POW(y-.5,2.0)), 0.5);
	    BaryonField[iBy][ELT(i,j,k)] = BField*(x-.5)/POW( (POW(x-.5,2.0)+POW(y-.5,2.0)), 0.5);
	    BaryonField[iBz][ELT(i,j,k)] = 0.0;
	  } else {
	    BaryonField[iBx][ELT(i,j,k)] = 
	      BaryonField[iBy][ELT(i,j,k)] = 
	      BaryonField[iBz][ELT(i,j,k)] = 0.0;

	  }

	} else if(FieldGeometry>1){
	  ENZO_FAIL("FieldGeometry > 1 not implemented yet!\n");
	}

      } // for(i...)  (loop over grid and set values)

  if (debug) {
    printf("Exiting ConductionTestInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;
}
