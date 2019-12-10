/***********************************************************************
/
/  GRID CLASS (RESTORE CONSISTENCY BETWEEN TOTAL AND INTERNAL ENERGY)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1: Tom Abel 2010 Added MHD
/
/  PURPOSE:
/    If using the dual energy formalism, this routine will restore
/     consistency between the total and internal (gas) energy fields.
/     This is done to either the entire field or just the boundary zones.
/     This can result in the supression of small shocks, so the routine
/     should be used with caution (i.e. just after interpolation which
/     generates errors of order a few percent anyway).
/
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
 
/* function prototypes */
 
 
/* InterpolateBoundaryFromParent function */
 
int grid::RestoreEnergyConsistency(int Region)
{
 
  int i, j, k, n;
 
  /* Error check */
 
  if (Region != ENTIRE_REGION && Region != ONLY_BOUNDARY) {
    ENZO_VFAIL("Region type %"ISYM" unknown.\n", Region)
  }
 
  /* If there is no work, we're done. */
 
  if (NumberOfBaryonFields <= 0 || DualEnergyFormalism == 0)
    return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num,B2Num,B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num,B2Num,B3Num) == FAIL) {
    fprintf(stderr, "\n");
    ENZO_FAIL("Error in grid->IdentifyPhysicalQuantities.");
  }
 
  /* a) Correct the entire field. */
 
  if (Region == ENTIRE_REGION) {
 
    /* Compute size of field. */
 
    int size = 1;
    for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Perform correction:  E = e + 1/2*v^2. */
 
    for (i = 0; i < size; i++)
      BaryonField[TENum][i] = BaryonField[GENum][i] +
	         0.5*(BaryonField[Vel1Num][i])*(BaryonField[Vel1Num][i]);
 
    if (GridRank > 1)
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] +=
	         0.5*(BaryonField[Vel2Num][i])*(BaryonField[Vel2Num][i]);
 
    if (GridRank > 2)
      for (i = 0; i < size; i++)
	BaryonField[TENum][i] +=
	         0.5*(BaryonField[Vel3Num][i])*(BaryonField[Vel3Num][i]);

    if (HydroMethod == MHD_RK) {
      float B2 = 0;
      for (i = 0; i < size; i++) {
	B2 = pow(BaryonField[B1Num][i],2) + pow(BaryonField[B2Num][i],2) +
	  pow(BaryonField[B3Num][i],2);
	BaryonField[TENum][i] += 0.5 * B2 / BaryonField[DensNum][i];
      }
    }
 
  } // end: Region == ENTIRE_FIELD
 
  /* b) Correct just the boundary zones. */
 
  if (Region == ONLY_BOUNDARY) {
 
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++)
	  if (i < GridStartIndex[0] || i > GridEndIndex[0] ||
	      j < GridStartIndex[1] || j > GridEndIndex[1] ||
	      k < GridStartIndex[2] || k > GridEndIndex[2]   ) {
 
	    n = k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i;
 
	    /* Perform correction:  E = e + 1/2*v^2. */
 
	    BaryonField[TENum][n] = BaryonField[GENum][n] +
	         0.5*(BaryonField[Vel1Num][n])*(BaryonField[Vel1Num][n]);
 
	    if (GridRank > 1)
	      BaryonField[TENum][n] +=
	         0.5*(BaryonField[Vel2Num][n])*(BaryonField[Vel2Num][n]);
 
	    if (GridRank > 2)
	      BaryonField[TENum][n] +=
	         0.5*(BaryonField[Vel3Num][n])*(BaryonField[Vel3Num][n]);
 
	    if (HydroMethod == MHD_RK) {

	      float B2; 
	      B2 = pow(BaryonField[B1Num][n],2) + pow(BaryonField[B2Num][n],2) +
		pow(BaryonField[B3Num][n],2);
	      BaryonField[TENum][n] += 0.5 * B2 / BaryonField[DensNum][n];
	    }


	  }
 
  } // end: Region == BOUNDARY_ONLY
 
  return SUCCESS;
 
}
