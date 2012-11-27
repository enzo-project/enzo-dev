/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS OUTFLOW BOUNDARY CONDITIONS FOR SHOCK POOL)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
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
#include "ShockPoolGlobalData.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
// Set the Left BoundaryValue of the chosen wave direction (set by
//  ShockPoolAngle) to the appropriate inflow boundary condition.
//  This is a sinusoidal wave
 
int ExternalBoundary::SetShockPoolBoundary(FLOAT time)
{
  /* declarations */
 
  int i, j, dim, index;
  int NumberOfZones[MAX_DIMENSION], Offset[MAX_DIMENSION];
  float deltime, distance, pos[MAX_DIMENSION];
  const float TwoPi = 6.283185;
 
  /* Compute size of entire mesh. */
 
  int size = 1;
  for (dim = 0; dim < BoundaryRank; dim++)
    size = size*BoundaryDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* set the appropriate BoundaryValues on the left side */
 
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
 
      /* If the BoundaryValue fields are missing, create them. */
 
      for (int field = 0; field < NumberOfBaryonFields; field++)
	if (BoundaryValue[field][dim][0] == NULL)
	  BoundaryValue[field][dim][0] =
	    new float[size/BoundaryDimension[dim]];
 
      /* Compute quantities needed for boundary face loop (below). */
 
      int dim1, dim2;
      dim1 = (dim == 0) ? 1 : 0;
      dim2 = dim1 + 1;
      dim2 = (dim2 == dim) ? dim2+1 : dim2;
      for (i = 0; i < 3; i++) {
	NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
	Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]) - 1;
      }
      pos[dim] = 0.0;
 
      /* Loop over the boundary face. */
 
      for (i = 0; i < BoundaryDimension[dim1]; i++)
	for (j = 0; j < BoundaryDimension[dim2]; j++) {
 
	  /* Compute the index into the boundary value. */
 
	  index = j*BoundaryDimension[dim1] + i;
 
	  /* Find the 3D vector from the corner to the current location. */
 
	  pos[dim1] = (float(i-Offset[dim1]))*
	    (DomainRightEdge[dim1]-DomainLeftEdge[dim1]) /
	      float(NumberOfZones[dim1]);
	  pos[dim2] = (float(j-Offset[dim2]))*
	    (DomainRightEdge[dim2]-DomainLeftEdge[dim2]) /
	      float(NumberOfZones[dim2]);
 
	  /* Compute the distance along the wave propogation vector
	     (cos(angle), sin(angle), 0). Convert to radians. */
	
	  distance = pos[0]*cos(ShockPoolAngle*TwoPi/360.0) +
	             pos[1]*sin(ShockPoolAngle*TwoPi/360.0);
 
	  /* Find the difference between the current time and the time at
	     which the wave will reach this point. */
 
	  deltime = time - distance/ShockPoolShockSpeed;
 
	  /* If deltime < 0, the shock has not yet reached this point. */
 
	  if (deltime > 0.0) {
 
	    /* Shock has arrived, set fields to post-shock values. */
 
	    BoundaryValue[DensNum][dim][0][index] = ShockPoolShockDensity;
	    BoundaryValue[TENum][dim][0][index] = ShockPoolShockTotalEnergy;
	    BoundaryValue[Vel1Num][dim][0][index] = ShockPoolShockVelocity[0];
	    if (BoundaryRank > 1)
	      BoundaryValue[Vel2Num][dim][0][index] =
	                                             ShockPoolShockVelocity[1];
	    if (BoundaryRank > 2)
	      BoundaryValue[Vel3Num][dim][0][index] =
	                                             ShockPoolShockVelocity[2];
 
	  } else {
 
	    /* If not, set the fields to their pre-shock values. */
 
	    BoundaryValue[DensNum][dim][0][index] = ShockPoolDensity;
	    BoundaryValue[TENum][dim][0]  [index] = ShockPoolTotalEnergy;
	    BoundaryValue[Vel1Num][dim][0][index] = ShockPoolVelocity[0];
	    if (BoundaryRank > 1)
	      BoundaryValue[Vel2Num][dim][0][index] = ShockPoolVelocity[1];
	    if (BoundaryRank > 2)

	      BoundaryValue[Vel3Num][dim][0][index] = ShockPoolVelocity[2];
 
	  }
 
	} // end loop over boundary slice
 
    } // end loop over boundary directions
 
  return SUCCESS;
 
}
