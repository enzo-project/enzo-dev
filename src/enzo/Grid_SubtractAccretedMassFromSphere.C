
/***********************************************************************
/
/  GRID: SUBTRACT ACCRETED MASS FROM CELLS
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
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
#include "Hierarchy.h"
#include "CosmologyParameters.h"

int FindField(int field, int farray[], int numfields);

int grid::SubtractAccretedMassFromSphere(Star *cstar, int level, float radius, float DensityUnits, 
					 float LengthUnits, float VelocityUnits, 
					 float TemperatureUnits, float TimeUnits, double EjectaDensity, 
					 int &CellsModified)
{

  const double Msun = 1.989e33;
  const double c = 3.0e10;

  int dim, i, j, k, index;
  FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];
  float coef, speed, maxVelocity;
  float OldDensity;
  float r1, norm, ramp, factor, newGE, increase;
  double old_mass;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (MBHAccretion <= 0) 
    return SUCCESS;

  /* Check if sphere overlaps with this grid */

  for (dim = 0; dim < GridRank; dim++)
    if (cstar->pos[dim] - radius > GridRightEdge[dim] ||
	cstar->pos[dim] + radius < GridLeftEdge[dim])
      return SUCCESS;

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Find Metallicity or SNColour field and set flag. */

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum; 
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    ENZO_FAIL("");
  }

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;


  /***********************************************************************
                            START SUBTRACTION
  ************************************************************************/

  float outerRadius2;

  // Correct if the volume with 27 cells is larger than the energy bubble volume
  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;
  if (BoxVolume > BubbleVolume) {
    //printf("Reducing ejecta density by %g\n", BubbleVolume / BoxVolume);
    EjectaDensity *= BubbleVolume / BoxVolume;
  }

  outerRadius2 = 1.2 * 1.2 * radius * radius;

  /* Update velocity of the star.  Because EjectaDensity will be subtracted out 
     with zero net momentum, increase the particle's velocity accordingly. 
     Note that the gas mass is already added to the star in Star_Accrete.C */

//  printf("grid::SAMFS: old_vel[1] = %g\n", cstar->vel[1]);
  old_mass = cstar->Mass + 
    EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  
  cstar->vel[0] *= old_mass / cstar->Mass; 
  cstar->vel[1] *= old_mass / cstar->Mass;
  cstar->vel[2] *= old_mass / cstar->Mass; 
//  printf("grid::SAMFS: old_mass = %lf  ->  cstar->Mass = %lf\n", old_mass, cstar->Mass);  
//  printf("grid::SAMFS: new_vel[1] = %g\n", cstar->vel[1]);

  for (k = 0; k < GridDimension[2]; k++) {
    
    delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
    delz = fabs(delz);
    delz = min(delz, DomainWidth[2]-delz);
    
    for (j = 0; j < GridDimension[1]; j++) {
      
      dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
      dely = fabs(dely);
      dely = min(dely, DomainWidth[1]-dely);

      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = 0; i < GridDimension[0]; i++, index++) {
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	delx = fabs(delx);
	delx = min(delx, DomainWidth[0]-delx);
	
	radius2 = delx*delx + dely*dely + delz*delz;
	if (radius2 <= outerRadius2) {
	  
	  r1 = sqrt(radius2) / radius;
	  norm = 0.98;
	  ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
	  //	     ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);
	  
	  /* 1/1.2^3 factor to dilute the density since we're
	     depositing a uniform ejecta in a sphere of 1.2*radius
	     without a ramp.  The ramp is only applied to the
	     energy*density factor. */
	  factor = 0.578704;
	  
	  OldDensity = BaryonField[DensNum][index];
	  BaryonField[DensNum][index] += factor*EjectaDensity;
	  increase = BaryonField[DensNum][index] / OldDensity;
//	  printf("grid::SAMFS: increase = %g\n", increase);

	  /* Update species and colour fields */
	  
	  if (MultiSpecies) {
	    BaryonField[DeNum][index] *= increase;
	    BaryonField[HINum][index] *= increase;
	    BaryonField[HIINum][index] *= increase;
	    BaryonField[HeINum][index] *= increase;
	    BaryonField[HeIINum][index] *= increase;
	    BaryonField[HeIIINum][index] *= increase;
	  }
	  if (MultiSpecies > 1) {
	    BaryonField[HMNum][index] *= increase;
	    BaryonField[H2INum][index] *= increase;
	    BaryonField[H2IINum][index] *= increase;
	  }
	  if (MultiSpecies > 2) {
	    BaryonField[DINum][index] *= increase;
	    BaryonField[DIINum][index] *= increase;
	    BaryonField[HIINum][index] *= increase;
	    BaryonField[HDINum][index] *= increase;
	  }
	  
	  if (MetallicityField == TRUE)
	    BaryonField[MetalNum][index] *= increase;
	  
	  if (MBHColourNum > 0)
	    BaryonField[MBHColourNum][index] *= increase;    
	  
	  CellsModified++;
	  
	} // END if inside radius
      }  // END i-direction
    }  // END j-direction
  }  // END k-direction
  
  return SUCCESS;

}
