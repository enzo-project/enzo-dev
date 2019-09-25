/***********************************************************************
/
/  GRID CLASS (INTERPOLATE PARTICLE POSITIONS FROM THE ACCELERATION FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: THIS ROUTINE DOES NOT TRANSFER DATA BETWEEN PROCESSORS!
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#define APDEBUG 0

/* function prototypes */
void ActiveParticleResetAccelerations(float *ActiveParticleAcceleration);
 
int grid::InterpolateParticlePositions(grid *FromGrid, int DifferenceType)
{
 
  FLOAT HoldLeftEdge[MAX_DIMENSION];
 
  /* Loop over all active dimensions */

  int dim, dim1;
 
  for (int dim = 0; dim < GridRank+ComputePotential; dim++) {
 
    /* Adjust the grid position if the acceleration is face-centered. */
 
    if (DifferenceType == DIFFERENCE_TYPE_STAGGERED &&
        dim != GridRank) {
      HoldLeftEdge[dim] = FromGrid->CellLeftEdge[dim][0];
      FromGrid->CellLeftEdge[dim][0] -= 0.5*FromGrid->CellWidth[dim][0];
    }

    if (NumberOfParticles > 0)
      if (FromGrid->InterpolatePositions(ParticlePosition, dim,
					 ParticleAcceleration[dim],
					 NumberOfParticles) == FAIL) {
        ENZO_FAIL("Error in grid->InterpolatePositions.\n");
      }

    if (NumberOfActiveParticles > 0) {

      FLOAT **ActiveParticlePosition = new FLOAT*[GridRank];
      for (dim1 = 0; dim1 < GridRank; dim1++)
        ActiveParticlePosition[dim1] = new FLOAT[NumberOfActiveParticles];
      this->GetActiveParticlePosition(ActiveParticlePosition);

      if (FromGrid->InterpolatePositions(ActiveParticlePosition, dim,
                                         ActiveParticleAcceleration[dim],
                                         NumberOfActiveParticles) == FAIL) {
        ENZO_FAIL("Error in grid->InterpolatePositions.\n");
      }
      /* Reset Active Particle positions if required */
      ActiveParticleResetAccelerations(ActiveParticleAcceleration[dim]);

      for (dim1 = 0; dim1 < GridRank; dim1++)
        delete [] ActiveParticlePosition[dim1];
      delete [] ActiveParticlePosition;
    }
 
      if(ProblemType==29){
        for(int i=0; i<NumberOfParticles; i++)
          printf("%"PISYM"  particle accelerations:  %e %e %e\n", ParticleNumber[i],
                 ParticleAcceleration[0][i],
                 ParticleAcceleration[1][i],
                 ParticleAcceleration[2][i]);
        fflush(stdout);
      }

      /* Adjust back. */
 
      if ( DifferenceType == DIFFERENCE_TYPE_STAGGERED && dim != GridRank)
        FromGrid->CellLeftEdge[dim][0] = HoldLeftEdge[dim];
  }
  
  return SUCCESS;
}
