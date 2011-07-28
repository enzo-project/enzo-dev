/***********************************************************************
/
/  GRID CLASS (COMPUTE PARTICLE AND GRID ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
 
/* function prototypes */
 
 
/* EvolveHierarchy function */
 
int grid::ComputeAccelerations(int level)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  int DiffType;

  if (NumberOfParticles > 0) {
  
    /* Set differencing type to be used (normal, face-centered) or
       staggered (cell-centered).  Staggered will generate a self-force. */

    DiffType = DIFFERENCE_TYPE_NORMAL;
    //    DiffType = DIFFERENCE_TYPE_STAGGERED;

    /* Compute the acceleration field for particles from potential. */
 
    if (this->ComputeAccelerationField(DiffType, level) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeAccelerationField.\n");
    }
 
    /* Add any fixed (external) acceleration to field. */
/*
    if (this->AddExternalAcceleration() == FAIL) {
      ENZO_FAIL("Error in grid->AddFixedAcceleration.\n");
    }
*/
    /* Clear particle accelerations. */
 
    this->DeleteParticleAcceleration();
    this->ClearParticleAccelerations();
 
    /* Move particles 1/2 step forward in preparation for interpolation. */
 
    this->UpdateParticlePosition(0.5*dtFixed);
 
    /* Interpolate the accelerations back to the grid and particles. */
 
    this->InterpolateParticlePositions(this, DiffType);
 
    /* Move particles 1/2 step backwards to return to original positions. */
 
    this->UpdateParticlePosition(-0.5*dtFixed);
 
    /* Clean up. */
 
    this->DeleteAccelerationField();
 
  } // end: if (NumberOfParticles > 0)
 
  /* Compute acceleration field for cells. */
 
  if (NumberOfBaryonFields > 0){

    /* Zeus hydro requires the acceleration field to be face-centered,
       while ppm assumed it is cell-centered. */

    DiffType = (HydroMethod == Zeus_Hydro) ? 
      DIFFERENCE_TYPE_STAGGERED : DIFFERENCE_TYPE_NORMAL;

     if (this->ComputeAccelerationField(DiffType, level) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeAccelerationField.\n");
    }

  } // end: if (NumberOfBaryonFields > 0)


  return SUCCESS;
}
