/***********************************************************************
/
/  GRID CLASS (RETURN PARTICLE INFO RELATIVE TO A GIVEN CENTER)
/
/  written by: Greg Bryan
/  date:       March, 2000
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/Fluxes.h"
#include "../enzo/GridList.h"
#include "../enzo/ExternalBoundary.h"
#include "../enzo/Grid.h"
#include "../enzo/CosmologyParameters.h"


/* function prototypes */



int grid::CollectParticleInfo(FLOAT SphereCenter[MAX_DIMENSION],
			      float SphereRadius,
			      int *ParticleCount, float *ParticleRadius,
			      float *ParticleDensity, float *ParticleVolume)
{  

  /* Return if the data is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int dim, i, j, k, index;
  float radius2, delx, dely, delz;
  FLOAT DomainWidth[MAX_DIMENSION];

  if (ComovingCoordinates != 1) {  
    InitialRedshift = 0; 
    //FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  } 

  /* Quick check to see if sphere overlaps this grid. */

  for (dim = 0; dim < GridRank; dim++)
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;

  /* Set domain width. */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  /* Find fields: density, total energy, velocity1-3. */

  if (NumberOfBaryonFields > 0) {

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);

    /* Loop over grid. */

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      if (GridRank > 1) {
	delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SphereCenter[2];
	delz = min(delz, DomainWidth[2]-delz);
      }
      else
	delz = 0;

      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	if (GridRank > 0) {
	  dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SphereCenter[1];
	  dely = min(dely, DomainWidth[1]-dely);
	}
	else
	  dely = 0;
	index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];

	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	  if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	    continue;
	
	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
	  delx = min(delx, DomainWidth[0]-delx);
	  
	  radius2 = delx*delx + dely*dely + delz*delz;

	  if (radius2 <= SphereRadius*SphereRadius) {
	    ParticleRadius[*ParticleCount] = sqrt(radius2);
	    ParticleDensity[*ParticleCount] = BaryonField[DensNum][index];
	    ParticleVolume[*ParticleCount] = pow(CellWidth[0][0], 3);
	    (*ParticleCount)++;
	  }

	} // end i loop
      } // end j loop
    } // end k loop

  } // end: if (NumberOfBaryonFields > 0)

  /* Loop over particles. */

  dely = delz = 0;
 
  for (i = 0; i < NumberOfParticles; i++) {

                      delx = ParticlePosition[0][i] - SphereCenter[0];
    if (GridRank > 1) dely = ParticlePosition[1][i] - SphereCenter[1];
    if (GridRank > 2) delz = ParticlePosition[2][i] - SphereCenter[2];

    radius2 = delx*delx + dely*dely + delz*delz;

    /* Add particle info to list (negative sign means particle). */

    if (radius2 <= SphereRadius*SphereRadius) {
      ParticleRadius[*ParticleCount] = sqrt(radius2);
      ParticleDensity[*ParticleCount] = -ParticleMass[i];
      ParticleVolume[*ParticleCount] = pow(CellWidth[0][0], 3);
      (*ParticleCount)++;
    }

  } // end: loop over particles

  return SUCCESS;
}
 
