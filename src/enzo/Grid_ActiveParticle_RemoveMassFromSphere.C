/***********************************************************************
/
/  GRID: REMOVE ACCRETED MASS FROM NEARBY CELLS
/
/  written by: Simone Gordon
/  date:       February, 2022
/  based on:   Grid_SubtractAccretedMassFromSphere.C
/              Grid_AddFeedbackSphere.C
/  modified1: 
/
/  PURPOSE: Subtract mass from gas grids within a pre-defined radius 
/           after forming ActiveParticle.
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
#include "phys_constants.h"

int FindField(int field, int farray[], int numfields);

int grid::RemoveMassFromSphere(ActiveParticleType* SS, 
           int level, float radius, float DensityUnits, 
					 float LengthUnits, float VelocityUnits, 
					 float TemperatureUnits, float TimeUnits, float Subtraction, 
					 int &CellsModified)
{

  int dim, i, j, k, index;
  FLOAT delx, dely, delz, radius2, DomainWidth[MAX_DIMENSION];
  double decrease;

  /* MY PROC = MPI PROCESS, Processor that owns this grid */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check if sphere overlaps with this grid */

  for (dim = 0; dim < GridRank; dim++)
    if (SS->pos[dim] - radius > GridRightEdge[dim] ||
	SS->pos[dim] + radius < GridLeftEdge[dim])
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
    MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
              MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;


  /***********************************************************************
                  Remove Mass From Grid After Formation
  ************************************************************************/

  /* Check if the volume with 3^3 cells is larger than the spherical volume 
     from which we remove the mass. If it is, scale down the Subtraction ratio. */

  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * pi / 3.0) * radius * radius * radius;
  if (BoxVolume > BubbleVolume) {
    fprintf(stderr, "%s: Level(%d) probably too coarse, rescaling Subtraction!\n",
                         __FUNCTION__, level);
    Subtraction *= BubbleVolume/BoxVolume;
  }

  /* Calculate how much the cell quantities are to be reduced by */
  decrease = max(1-Subtraction, 0.5);

  for (k = 0; k < GridDimension[2]; k++) {
    
    delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SS->pos[2];
    delz = fabs(delz);
    delz = min(delz, DomainWidth[2]-delz);
    
    for (j = 0; j < GridDimension[1]; j++) {
      
      dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SS->pos[1];
      dely = fabs(dely);
      dely = min(dely, DomainWidth[1]-dely);

      index = (k*GridDimension[1] + j)*GridDimension[0];
      for (i = 0; i < GridDimension[0]; i++, index++) {

        delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SS->pos[0];
        delx = fabs(delx);
        delx = min(delx, DomainWidth[0]-delx);
        
        radius2 = delx*delx + dely*dely + delz*delz;
        if (radius2 <= radius*radius) {

          /* Update density */

          float density1 = BaryonField[DensNum][index];

          BaryonField[DensNum][index] *= decrease;

          float density2 = BaryonField[DensNum][index];

          /* Update velocities and TE. The grid lost some mass, so velocity is increased.
             For DualEnergyFormalism = 0, you don't have to update any energy field */

          if (GENum >= 0 && DualEnergyFormalism) 
            for (dim = 0; dim < GridRank; dim++)
              BaryonField[TENum][index] -= 
              0.5 * BaryonField[Vel1Num+dim][index] * 
              BaryonField[Vel1Num+dim][index];
          
              BaryonField[Vel1Num][index] /= decrease;
              BaryonField[Vel2Num][index] /= decrease;
              BaryonField[Vel3Num][index] /= decrease;
          
          if (GENum >= 0 && DualEnergyFormalism) 
            for (dim = 0; dim < GridRank; dim++)
              BaryonField[TENum][index] += 
              0.5 * BaryonField[Vel1Num+dim][index] * 
              BaryonField[Vel1Num+dim][index];
          
          /* Update species and colour fields */
          
          if (MultiSpecies) {
            BaryonField[DeNum][index] *= decrease;
            BaryonField[HINum][index] *= decrease;
            BaryonField[HIINum][index] *= decrease;
            BaryonField[HeINum][index] *= decrease;
            BaryonField[HeIINum][index] *= decrease;
            BaryonField[HeIIINum][index] *= decrease;
          }
          if (MultiSpecies > 1) {
            BaryonField[HMNum][index] *= decrease;
            BaryonField[H2INum][index] *= decrease;
            BaryonField[H2IINum][index] *= decrease;
          }
          if (MultiSpecies > 2) {
            BaryonField[DINum][index] *= decrease;
            BaryonField[DIINum][index] *= decrease;
            BaryonField[HDINum][index] *= decrease;
          }
          
          if (MetallicityField == TRUE)
            BaryonField[MetalNum][index] *= decrease;
          
          
          /* Increment number of cells modified */
          CellsModified++;
          
        } // END if inside radius
      }  // END i-direction
    }  // END j-direction
  }  // END k-direction


  fprintf(stderr, "%s: Mass removed from radius = %e pc, cells modified = %"ISYM", fractional decrease = %e.\n", 
          __FUNCTION__,
          radius * LengthUnits / pc_cm,
          CellsModified,
          decrease);
  
  return SUCCESS;

}
