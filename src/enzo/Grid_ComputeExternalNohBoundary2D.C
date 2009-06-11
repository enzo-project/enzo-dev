/***********************************************************************
/
/  GRID CLASS (COMPUTE RIGHT BOUNDARIES FOR 2D NOH PROBLEM BASED ON 
/              THE EXACT ANALYTICAL SOLUTION)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, May 2005 - Nov 2005; Feb 2007
/
/  PURPOSE:
/    This function computes boundary values for the external boundaries
/    using analytical formulae:
/        density = (1 + t/sqrt(x^2 + y^2))          2D-case
/        pressure = 1e-6
/        velocity = 1 (radial)
/
/    2D test is the only only supported.
/
/  RETURNS: FAIL or SUCCESS
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

/* ComputeExternalNohBoundary function. It is assumed that subgrids do not 
   touch the time-dependent external boundary, so this function will not
   operate on AMR subgrids. */

int grid::ComputeExternalNohBoundary()
{

  float d0 = 1.0, p0 = 1.0e-6, u0 = -1.0;

  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
    //    ERROR_MESSAGE; // THIS ROUTINE MUST(?) BE MODIFIED TO SUPPORT DOMAIN DECOMPOSITION

  if (NumberOfBaryonFields > 0) {


  /* Compute offsets from right and left domain edges. */

  int dim, GridOffsetLeft[MAX_DIMENSION], GridOffsetRight[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (dim < GridRank) {
      GridOffsetLeft[dim]  = nint(( GridLeftEdge[dim] - DomainLeftEdge[dim])/
                             CellWidth[dim][0]);
      GridOffsetRight[dim] = nint((GridRightEdge[dim] - DomainRightEdge[dim])/
                             CellWidth[dim][0]);
    }
    else {
      GridOffsetLeft[dim]  = 0;
      GridOffsetRight[dim] = 0;
    }

    /* Compute needed portion of the field on current grid,
       only if this grid is in contact with the domain boundary. 
       Do NOT include (0,1) and (1,0)-corners of the domain (if
       only the upper right quadrant is used, not the full box), but
       include corners for those subgrids which do not cover any 
       of those domain corners.
    */
    
    int i, j, index, 
      ishift = GridStartIndex[0], 
      jshift = GridStartIndex[1];

    if (NohProblemFullBox == 1) {
      if ( GridDimension[0]%2 != 0 || 
	  (GridDimension[1]%2 != 0 && GridRank > 1))
	ERROR_MESSAGE; // this must be eventually redone; works for serial runs only
      ishift += nint((DomainRightEdge[0] - DomainLeftEdge[0])/CellWidth[0][0])/2;
      jshift += nint((DomainRightEdge[1] - DomainLeftEdge[1])/CellWidth[1][0])/2;
    }

     int istart = 0, jstart = 0;
    if (GridOffsetLeft[0] == 0 && NohProblemFullBox == 0)
      istart = GridStartIndex[0];
    else ishift -= GridOffsetLeft[0];

    if (GridOffsetLeft[1] == 0 && NohProblemFullBox == 0)
      jstart = GridStartIndex[1];
    else jshift -= GridOffsetLeft[1];

    float radius, xx, yy;
    float time = Time;

    /* j-sweeps, LEFT boundary. */

    if (GridOffsetLeft[0] == 0 && NohProblemFullBox == 1)
      for (j = jstart; j < GridDimension[1]; j++) {
	index = j*GridDimension[0];
	yy = j + 0.5 - jshift;
	for (i = 0; i < GridStartIndex[0]; i++) {
	  xx = i + 0.5 - ishift;
	  radius = max(tiny_number, sqrt(xx*xx + yy*yy));
	  BaryonField[0][index+i] = d0 + time/radius/CellWidth[0][0];
	  BaryonField[1][index+i] = p0/(Gamma-1.0)/d0;
	  BaryonField[2][index+i] = u0*xx/radius;
	  BaryonField[3][index+i] = u0*yy/radius;
	  if (HydroMethod != Zeus_Hydro)
	    BaryonField[1][index+i] +=  0.5*u0*u0;
	}
      }

    /* i-sweeps, LEFT boundary. */

    if (GridOffsetLeft[1] == 0 && NohProblemFullBox == 1)
      for (j = 0; j < GridStartIndex[1]; j++) {
	index = j*GridDimension[0];
	yy = j + 0.5 - jshift;
	for (i = istart; i < GridDimension[0]; i++) {
	  xx = i + 0.5 - ishift;
	  radius = max(tiny_number, sqrt(xx*xx + yy*yy));
	  BaryonField[0][index+i] = d0 + time/radius/CellWidth[0][0];
	  BaryonField[1][index+i] = p0/(Gamma-1.0)/d0;
	  BaryonField[2][index+i] = u0*xx/radius;
	  BaryonField[3][index+i] = u0*yy/radius;
	  if (HydroMethod != Zeus_Hydro)
	    BaryonField[1][index+i] += 0.5*u0*u0;
	}
      }      

    /* j-sweeps, only RIGHT boundary. */

    if (GridOffsetRight[0] == 0)
      for (j = jstart; j < GridDimension[1]; j++) {
	index = j*GridDimension[0];
	yy = j + 0.5 - jshift;
	for (i = GridEndIndex[0]+1; i < GridDimension[0]; i++) {
	  xx = i + 0.5 - ishift;
	  radius = max(tiny_number, sqrt(xx*xx + yy*yy));
	  BaryonField[0][index+i] = d0 + time/radius/CellWidth[0][0];
	  BaryonField[1][index+i] = p0/(Gamma-1.0)/d0;
	  BaryonField[2][index+i] = u0*xx/radius;
	  BaryonField[3][index+i] = u0*yy/radius;
	  if (HydroMethod != Zeus_Hydro)
	    BaryonField[1][index+i] +=  0.5*u0*u0;
	}
      }
      

    /* i-sweeps, only RIGHT boundary. */

    if (GridOffsetRight[1] == 0)
      for (j = GridEndIndex[1]+1; j < GridDimension[1]; j++) {
	index = j*GridDimension[0];
	yy = j + 0.5 - jshift;
	for (i = istart; i < GridDimension[0]; i++) {
	  xx = i + 0.5 - ishift;
	  radius = max(tiny_number, sqrt(xx*xx + yy*yy));
	  BaryonField[0][index+i] = d0 + time/radius/CellWidth[0][0];
	  BaryonField[1][index+i] = p0/(Gamma-1.0)/d0;
	  BaryonField[2][index+i] = u0*xx/radius;
	  BaryonField[3][index+i] = u0*yy/radius;
	  if (HydroMethod != Zeus_Hydro)
	    BaryonField[1][index+i] += 0.5*u0*u0;
	}
      }

    /* DualEnergyFormalism not supported. */

    if (DualEnergyFormalism)
	ERROR_MESSAGE;

  } // end: if (NumberOfBaryonFields > 0)

  this->DebugCheck("ComputeExternalNohBoundary (after)");
  
  return SUCCESS;

}
