/***********************************************************************
/
/  Apply gravity to the temporary grid for Galaxy Particles.
/
/  written by: Stephen Skory
/  date:       October, 2012
/
/  note:       
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "ActiveParticle_GalaxyParticle.h"

#define NO_DEBUG

int grid::ApplyGalaxyParticleGravity(ActiveParticleType** ThisParticle){

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;

  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

   FLOAT *ParticlePosition = (*ThisParticle)->ReturnPosition();
   float rad = static_cast<ActiveParticleType_GalaxyParticle*>(*ThisParticle)->Radius;
   
   if ((GridLeftEdge[0] > ParticlePosition[0]+rad) || (GridRightEdge[0] < ParticlePosition[0]-rad) ||
       (GridLeftEdge[1] > ParticlePosition[1]+rad) || (GridRightEdge[1] < ParticlePosition[1]-rad) ||
       (GridLeftEdge[2] > ParticlePosition[2]+rad) || (GridRightEdge[2] < ParticlePosition[2]-rad))
     return SUCCESS;

  
  int i, j, k, index;
  float rad2dx, dist2;
  float dx = float(this->CellWidth[0][0]);
  float max_gmf = 0.;
  FLOAT xx, yy, zz;

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum, MetalIINum, MetalAGBNum, MetalNSMNum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum,
                                 MetalIINum, MetalAGBNum, MetalNSMNum, MBHColourNum, 
                                 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");
  }
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
         Vel3Num, TENum) == FAIL) {
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
   }

  rad2dx = (rad + dx) * (rad + dx);
  
  printf("I would be applying gravity to this grid:\n");
  printf("%f %f %f : %f %f %f\n", GravitatingMassFieldLeftEdge[0],
    GravitatingMassFieldLeftEdge[1], GravitatingMassFieldLeftEdge[2],
    GravitatingMassFieldLeftEdge[0] + GravitatingMassFieldDimension[0]* GravitatingMassFieldCellSize,
    GravitatingMassFieldLeftEdge[1] + GravitatingMassFieldDimension[1]* GravitatingMassFieldCellSize,
    GravitatingMassFieldLeftEdge[2] + GravitatingMassFieldDimension[2]* GravitatingMassFieldCellSize);
//   printf("%d %d %d : %d %d %d\n", this->GridStartIndex[0], this->GridStartIndex[1],
//     this->GridStartIndex[2], this->GridEndIndex[0], this->GridEndIndex[1],
//     this->GridEndIndex[2]);
//   printf("%f %f %f\n", this->CellLeftEdge[0][3] + 0.5*this->CellWidth[0][3],
//     this->CellLeftEdge[1][3] + 0.5*this->CellWidth[1][3],
//     this->CellLeftEdge[2][3] + 0.5*this->CellWidth[2][3]);
  
    for (k = 0; k < GravitatingMassFieldDimension[0]; k++) {
      for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
        index = GRIDINDEX_NOGHOST(0, j, k);
        for (i = 0; i < this->GravitatingMassFieldDimension[0]; i++, index++) {
            max_gmf = max(this->GravitatingMassField[index], max_gmf);
        } // i
      } // j
    } // k
    // We actually need the distance from a cell to the particle, so these
    // heavy loops are required.
    for (k = 0; k < GravitatingMassFieldDimension[0]; k++) {
      zz = this->GravitatingMassFieldLeftEdge[2] + (float(k) + 0.5)*this->CellWidth[2][k];
      for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
        yy = this->GravitatingMassFieldLeftEdge[1] + (float(j) + 0.5)*this->CellWidth[1][j];
        index = GRIDINDEX_NOGHOST(0, j, k);
        for (i = 0; i < this->GravitatingMassFieldDimension[0]; i++, index++) {
          xx = this->GravitatingMassFieldLeftEdge[0] + (float(i) + 0.5)*this->CellWidth[0][i];
          // Don't need to consider periodicity because this grid is contiguous.
          dist2 = ((xx - ParticlePosition[0]) * (xx - ParticlePosition[0]) + 
                   (yy - ParticlePosition[1]) * (yy - ParticlePosition[1]) + 
                   (zz - ParticlePosition[2]) * (zz - ParticlePosition[2]));
          if (dist2 > rad2dx) continue;
          // If we've gotten this far, at least some part of this cell overlaps
          // with the sphere of this GP.
          this->GravitatingMassField[index] = 5 * max_gmf;
        } // i
      } // j
    } // k

  return SUCCESS;
}

#undef DEBUG
