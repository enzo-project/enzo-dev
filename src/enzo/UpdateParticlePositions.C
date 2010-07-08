/***********************************************************************
/
/  UPDATE PARTICLE POSITIONS
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
 
int UpdateParticlePositions(grid *Grid)
{

  float dt = Grid->ReturnTimeStep();

  LCAPERF_START("UpdateParticlePositions");

  /* 1) v(n) --> v(n+1/2) with a(n+1/2) */
 
  Grid->DebugCheck("UpdateParticlePosition step 1");
  if (Grid->UpdateParticleVelocity(0.5*dt) == FAIL) {
    ENZO_FAIL("Error in grid->UpdateParticleVelocity./\n");
  }
 
  /* If there are any tracer particles (which follow fluid flow),
     then set their velocities now (i.e. just before particle push). */
 
  if (TracerParticleOn)
    Grid->TracerParticleSetVelocity();
 
  /* 2) x(n) --> x(n+1) with v(n+1/2)
        (problem 23 is TestGravity, so don't move the particles). */
 
  Grid->DebugCheck("UpdateParticlePosition step 2");
  if (ProblemType != 23)
    if (Grid->UpdateParticlePosition(dt) == FAIL) {
      ENZO_FAIL("Error in grid->UpdateParticlePosition./\n");
    }
 
  /* 3) v(n+1/2) --> v(n+1) with a(n+1/2) */
 
  Grid->DebugCheck("UpdateParticlePosition step 3");
  if (Grid->UpdateParticleVelocity(0.5*dt) == FAIL) {
    ENZO_FAIL("Error in grid->UpdateParticleVelocity./\n");

  }
 
 
  LCAPERF_STOP("UpdateParticlePositions");
  return SUCCESS;
}
