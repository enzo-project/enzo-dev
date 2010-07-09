/***********************************************************************
/
/  GRID CLASS (ADD FIXED ACCELERATION TO ACCEL FIELDS)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE: Certain problems required external acceleration fields.
/    This routine adds them to the existing self-gravitating fields.
/
/  NOTE: This routine has really been replaced by
/        grid::ComputeAccelerationFieldExternal.  It should be removed.
/
************************************************************************/
 
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
#include "SphericalInfall.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
int grid::AddExternalAcceleration()
{
 
  /* -----------------------------------------------------------------
     Problem: SphericalInfall
     ----------------------------------------------------------------- */
 
  if (ProblemType == 24) {
 
    float Pi = 3.14159, accel;
    FLOAT a = 1, dadt, rcubed, xpos, ypos = 0, zpos = 0;
    int   i, j, k, n = 0;
 
    if (SphericalInfallFixedAcceleration) {
 
      /* Compute adot/a at time = t+1/2dt (time-centered). */
 
      if (ComovingCoordinates)
	if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	    == FAIL) {
	  ENZO_FAIL("Error in CosmologyComputeExpansionFactor.\n");
	}
 
      /* Loop over grid, adding acceleration to field. */
 
      for (k = 0; k < GridDimension[2]; k++) {
	if (GridRank > 2)
	  zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] -
	         SphericalInfallCenter[2];
	for (j = 0; j < GridDimension[1]; j++) {
	  if (GridRank > 1)
	    ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -
	           SphericalInfallCenter[1];
	  for (i = 0; i < GridDimension[0]; i++, n++) {
	    xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] -
	           SphericalInfallCenter[0];
 
	    /* Compute distance from center. */
 
	    rcubed = fabs(xpos*xpos*xpos) + fabs(ypos*ypos*ypos) +
	             fabs(zpos*zpos*zpos);
            if (rcubed < 1.0e-12)
	       printf("r3 = %"GOUTSYM" pos:%"GOUTSYM" %"GOUTSYM" %"GOUTSYM" Mid:%"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n", rcubed,
                      xpos, ypos, zpos, SphericalInfallCenter[0],
                      SphericalInfallCenter[1], SphericalInfallCenter[2]);
	    rcubed = max(rcubed, 1.0e-9);
 
	    /* Compute force. */
	    /* Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm().
	       (i.e. 1/a^2 * a = 1/a). */
 
	    accel = GravitationalConstant*SphericalInfallFixedMass/
	      (4.0*Pi*rcubed*a);
 
	    /* Apply force. */
 
	    AccelerationField[0][n] -= accel*xpos;
	    if (GridRank > 1)
	      AccelerationField[1][n] -= accel*ypos;
	    if (GridRank > 2)
	      AccelerationField[2][n] -= accel*zpos;
 
	  }
	}
      } // end: loop over grid
 
      /* DO PARTICLES! */
 
    } // end: if (SphericalInfallFixedAcceleration)

 
  } // end: SphericalInfall
 
 
  return SUCCESS;
}
 
