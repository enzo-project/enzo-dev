/***********************************************************************
/
/  GRID CLASS (ADD FIXED ACCELERATION TO ACCEL FIELDS)
/
/  written by: Greg Bryan
/  date:       June, 1996
/  modified1:
/
/  PURPOSE: Certain problems required external acceleration fields.
/    This routine adds them to the existing self-gravitating fields, or
/     creates new fields if necessary.
/    There is currently support for:
/      1) A uniform gravitational field in one of the orthogonal directions
/         for grid cells only
/      2) A 3D point source field for grid cells only
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
int grid::ComputeAccelerationFieldExternal()
{
 
  /* Return if this does not concern us */
  if (!(UniformGravity || PointSourceGravity || ExternalGravity)) return SUCCESS;

  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  int dim, i, j, k, size = 1;
 
  JBPERF_START("grid_ComputeAccelerationFieldExternal");

  /* Compute field size (in floats). */
 
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Check if acceleration field exists.  If not create it and zero it. */
 
  if (AccelerationField[0] == NULL)
    for (dim = 0; dim < GridRank; dim++) {
      AccelerationField[dim] = new float[size];
      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = 0;
    }
 
  /* -----------------------------------------------------------------
     Point Source gravity
     ----------------------------------------------------------------- */
 
  if (PointSourceGravity) {
 
    FLOAT a = 1, accel, dadt, rcubed, xpos, ypos = 0, zpos = 0, rcore,x ;
 
    /* Compute adot/a at time = t+1/2dt (time-centered). */
 
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	  == FAIL) {
	fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
	ENZO_FAIL("");
      }
 
    /* Loop over grid, adding acceleration to field. */
 
    for (dim = 0; dim < GridRank; dim++) {
      int n = 0;
 
      for (k = 0; k < GridDimension[2]; k++) {
	if (GridRank > 2)
	  zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] -
	    PointSourceGravityPosition[2];
	if (dim == 2 && HydroMethod == Zeus_Hydro)
	  zpos -= 0.5*CellWidth[2][k];
 
	for (j = 0; j < GridDimension[1]; j++) {
	  if (GridRank > 1)
	    ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -
	      PointSourceGravityPosition[1];
	  if (dim == 1 && HydroMethod == Zeus_Hydro)
	    ypos -= 0.5*CellWidth[1][j];
 
	  for (i = 0; i < GridDimension[0]; i++, n++) {
	    xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] -
	      PointSourceGravityPosition[0];
	    if (dim == 0 && HydroMethod == Zeus_Hydro)
	      xpos -= 0.5*CellWidth[0][i];
	
	    /* Compute distance from center. */
 
	    // rcubed = fabs(xpos*xpos*xpos) + fabs(ypos*ypos*ypos) +
	    //   fabs(zpos*zpos*zpos);
	    //
	    // expression below changed 9th November 2004  (BWO/RH)
 
	    rcubed = POW( (xpos*xpos + ypos*ypos  + zpos*zpos), 1.5);
 
	    rcore = max(0.1*CellWidth[0][0], PointSourceGravityCoreRadius);
 
	    /* (1) is a real (softened) point-source, (2) is NFW profile */
 
	    if (PointSourceGravity == 1)
	      rcubed += rcore*rcore*rcore;
	    else {
	      x = POW(rcubed, FLOAT(1.0/3.0))/rcore;
	      rcubed /= (log(1+x) - x/(1+x));
	    }
	
	    /* Compute force. */
	    /* Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm().
	       (i.e. 1/a^2 * a = 1/a). */
	
	    accel = PointSourceGravityConstant/(rcubed*a);
	
	    /* Apply force. */
	
	    if (dim == 0)
	      AccelerationField[0][n] -= accel*xpos;
	    if (dim == 1)
	      AccelerationField[1][n] -= accel*ypos;
	    if (dim == 2)
	      AccelerationField[2][n] -= accel*zpos;
 
	  }
	}
      } // end: loop over grid
    } // end: loop over dims
 
    /* DO PARTICLES HERE! */
 
  } // end: if (PointSourceGravity)
 
  /* -----------------------------------------------------------------
     Uniform gravity field
     ----------------------------------------------------------------- */
 
  if (UniformGravity) {
 
//    if (debug)
//      printf("grid->ComputeAccelerationFieldExternal: dir, g = %"ISYM", %"GSYM"\n",
//	     UniformGravityDirection, UniformGravityConstant);
		
    for (dim = 0; dim < GridRank; dim++) {
 
      /* Set constant for this dimension. */
 
      float Constant = 0.0;
      if (dim == UniformGravityDirection)
	Constant = UniformGravityConstant;
	
      /* Set field. */
 
      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = Constant;
 
    } // loop over dims
 
  } // end: if (UniformGravity)
 
  JBPERF_STOP("grid_ComputeAccelerationFieldExternal");
  return SUCCESS;
}
 
