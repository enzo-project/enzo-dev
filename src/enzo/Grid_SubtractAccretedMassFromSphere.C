/***********************************************************************
/
/  GRID: SUBTRACT ACCRETED MASS FROM NEARBY CELLS
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/  PURPOSE: Subtract mass from gas grids after accreting onto MBH,
/           For BlackHole, a different approach is used at the moment;
/           see Star_SubtractAccretedMassFromCell.C
/ 
/           This file is based on the logic of Grid_AddFeedbackSphere.C
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
					 float TemperatureUnits, float TimeUnits, double Subtraction, 
					 int &CellsModified)
{

  const double Msun = 1.989e33;
  const double c = 3.0e10;

  int dim, i, j, k, index;
  FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];
  float coef, speed, maxVelocity;
  float OldDensity;
  float r1, norm, ramp, factor, newGE;
  double old_mass, increase;

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

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;


  /***********************************************************************
                  MASS SUBTRACTION AFTER ACCRETION ONTO MBH
  ************************************************************************/

  // Correct if the volume with 27 cells is larger than the bubble volume 
  // from which we subtract the mass
  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;
  if (BoxVolume > BubbleVolume) {
    fprintf(stdout, "grid::SAMFS: level(%d) probably too coarse, rescaling Subtraction!\n", level);
    Subtraction *= BubbleVolume/BoxVolume;
  }

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
	if (radius2 <= radius*radius) {

	  increase = max(1-Subtraction, 0.9); 

#ifdef SUBTRACTION_UNIFORM 
	  // check CalculateSubtractionParameters.C if you want this
	  increase = max(BaryonField[DensNum][index] + Subtraction, 0.90*BaryonField[DensNum][index]) 
	    / BaryonField[DensNum][index];
#endif

	  /* Update density */
	  BaryonField[DensNum][index] *= increase;
	  // this "increase" method could be potentially problematic when the accretion rate is too low
//	  printf("grid::SAMFS: increase = %lf\n", increase); 

	  /* Update velocities and TE; the grid lost some mass, increasing the velocity;
	     for DualEnergyFormalism = 0 you don't have to update any energy field */

	  if (GENum >= 0 && DualEnergyFormalism) 
	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[TENum][index] -= 
		0.5 * BaryonField[Vel1Num+dim][index] * 
		BaryonField[Vel1Num+dim][index];
	  
	  BaryonField[Vel1Num][index] /= increase;
	  BaryonField[Vel2Num][index] /= increase;
	  BaryonField[Vel3Num][index] /= increase;
	  
	  if (GENum >= 0 && DualEnergyFormalism) 
	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[TENum][index] += 
		0.5 * BaryonField[Vel1Num+dim][index] * 
		BaryonField[Vel1Num+dim][index];
	  
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


//  printf("grid::SAMFS: radius (pc) = %lf, increase = %lf, mass subtracted (Msun) = %lf\n", 
//	 radius * LengthUnits / 3.086e18, increase, 
//	 Subtraction * (4*M_PI/3.0 * pow(radius*LengthUnits, 3)) * DensityUnits / Msun); 
  
  return SUCCESS;

}
