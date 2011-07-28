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
int grid::ConductionCloudInitialize (float CloudOverdensity, FLOAT CloudWidth, int CloudType) {

  if (debug) {
    printf("Entering ConductionCloudInitialize\n");
    fflush(stdout);
  }

  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;
  
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  FLOAT sig2 = CloudWidth*CloudWidth;

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
	r2 = POW(x-0.5, 2.0);

	if(GridRank>1){
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  r2 += POW(y-0.5, 2.0);
	}

	if(GridRank>2){
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  r2 += POW(z-0.5, 2.0);
	}

	celldist = POW(r2,0.5);

	float val=1.0;

	/* now we select between our different pulse options */

	if(CloudType == 1){  // gaussian pulse

	  val = 1.0 + exp(-1.0*r2/sig2/2.0)*(CloudOverdensity-1.0);

	} else if(CloudType == 2){ // square pulse

	  if(r2 <= sig2)
	    val = CloudOverdensity;
	  
	} else if(CloudType == 3){  // sinusoidal pulse with values along x-axis

	  val = 1.0 + CloudOverdensity + CloudOverdensity * sin(2.0 * 3.14158 * x / CloudWidth);
	  
	} else if(CloudType == 4){  // square pulse with smoothed edges (as suggested by A. Kravtsov)

	  val = 1.0 + (CloudOverdensity-1.0)*(1.0 - tanh((10.*(celldist/CloudWidth-1.0)))) / 2.0;

	} else {

	  ENZO_FAIL("Grid::ConductionCloudInitialize: CloudType is not 1,2, 3 or 4!");
	  
	}

	/*--------- MODIFY THE DENSITY HERE -----------*/
	BaryonField[DensNum][ELT(i,j,k)] /= val;

	if(HydroMethod==Zeus_Hydro){  // ZEUS
	  BaryonField[TENum][ELT(i,j,k)] *= val;  // TE = gas energy
	} else{ // PPM
	  
	  BaryonField[TENum][ELT(i,j,k)] *= val;  // TE = total energy energy, but velocity=0 here.

	  if(DualEnergyFormalism)
	    BaryonField[GENum][ELT(i,j,k)] *= val;  // if DEF=1, need to separately set the gas internal energy.
	}

	if(TestProblemData.MultiSpecies>1){
	  fprintf(stderr,"This problem type is not set up for MultiSpecies > 1.  Oops!\n");
	  ENZO_FAIL("Error in Grid::ConductionCloudInitialize.");
	}

	// Set multispecies fields!
	// this attempts to set them such that species conservation is maintained,
	// using the method in CosmologySimulationInitializeGrid.C
	if(TestProblemData.MultiSpecies) {

	  BaryonField[HIINum][ELT(i,j,k)] = TestProblemData.HII_Fraction * 
	    BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.HydrogenFractionByMass;
	  
	  BaryonField[HeIINum][ELT(i,j,k)] = TestProblemData.HeII_Fraction *
	    BaryonField[DensNum][ELT(i,j,k)]*(1.-TestProblemData.HydrogenFractionByMass);
	      
	  BaryonField[HeIIINum][ELT(i,j,k)] = TestProblemData.HeIII_Fraction *
	    BaryonField[DensNum][ELT(i,j,k)]*(1.-TestProblemData.HydrogenFractionByMass);

	  BaryonField[HeINum][ELT(i,j,k)] = BaryonField[DensNum][ELT(i,j,k)]*(1.-TestProblemData.HydrogenFractionByMass) -
	    BaryonField[HeIINum][ELT(i,j,k)] - BaryonField[HeIIINum][ELT(i,j,k)];

	  // HI density is calculated by subtracting off the various ionized fractions
	  // from the total
	  BaryonField[HINum][ELT(i,j,k)] = TestProblemData.HydrogenFractionByMass*BaryonField[DensNum][ELT(i,j,k)]
	    - BaryonField[HIINum][ELT(i,j,k)];
	  
	  // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
	  // density for convenience) is calculated by summing up all of the ionized species.
	  // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
	  // calculating mass density, not number density (because the BaryonField values are 4x as
	  // heavy for helium for a single electron)
	  BaryonField[DeNum][ELT(i,j,k)] = BaryonField[HIINum][ELT(i,j,k)] +
	    0.25*BaryonField[HeIINum][ELT(i,j,k)] + 0.5*BaryonField[HeIIINum][ELT(i,j,k)];
	  
	} // if(TestProblemData.MultiSpecies)

	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
	  BaryonField[MetalNum][ELT(i,j,k)] = 
	    BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;

      } // for(i...)  (loop over grid and set values)

  if (debug) {
    printf("Exiting ConductionCloudInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;
}
