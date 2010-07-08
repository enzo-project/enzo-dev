/***********************************************************************
/
/  GRID CLASS (RE-INITIALIZE FOR A SUPERNOVA EXPLOSION RUN)
/
/  written by: Greg Bryan
/  date:       February, 2000
/  modified1:
/
/  PURPOSE:
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
 
int grid::SupernovaRestartInitialize(float EjectaDensity, float EjectaRadius,
				     float EjectaThermalEnergy,
				     FLOAT EjectaCenter[3], int ColourField,
                                     int *NumberOfCellsSet)
{
  /* declarations */
 
  int dim, i, j, k, n;
  FLOAT delx, dely, delz, radius2, DomainWidth[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
 
  /* Set up colour field. */
 
  if (ColourField) {
    ENZO_FAIL("ColourField not implemented yet.\n");
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Loop over grid and set quantities. */
 
  n = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    if (GridRank > 1) {
      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - EjectaCenter[2];
      delz = min(delz, DomainWidth[2]-delz);
    }
    else
      delz = 0;
 
    for (j = 0; j < GridDimension[1]; j++) {
      if (GridRank > 0) {
	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - EjectaCenter[1];
	dely = min(dely, DomainWidth[1]-dely);
      }
      else
	dely = 0;
 
      for (i = 0; i < GridDimension[0]; i++, n++) {
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - EjectaCenter[0];
	delx = min(delx, DomainWidth[0]-delx);
 
	/* Compute square of distance from cell to center. */
 
	radius2 = delx*delx + dely*dely + delz*delz;
 
	if (radius2 <= EjectaRadius*EjectaRadius*1.2*1.2) {
 
	  float r1 = sqrt(radius2)/EjectaRadius;
	  float ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);
 
	  /* Cell is within ejecta, so set density, etc. */
 
	  BaryonField[DensNum][n] = EjectaDensity;
	  BaryonField[TENum][n]   = EjectaThermalEnergy*ramp;
	  if (GENum >= 0)
	    BaryonField[GENum][n]   = BaryonField[TENum][n];
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[Vel1Num+dim][n] = 0;
	  (*NumberOfCellsSet)++;
 
	} // end: if (radius2 < EjectaRadius*EjectaRadius)

 
      } // next i
    } // next j
  } // next j
 
 
  return SUCCESS;
}
