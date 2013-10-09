/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE PRESENCE OF SHOCKS)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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
 
int grid::FlagCellsToBeRefinedByShocks()
{
  /* declarations */
 
  int i, j, k, index, dim;
  float DelPressure, DelVelocity;
  float EnergyThisCell, EnergyRatio, Energy1, Energy2, Energy3;
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the pressure. */
 
  float *Pressure = new float[size];
  this->ComputePressure(Time, Pressure);
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* loop over active dimensions */
 
  int Offset = 1;
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* For shock: (p(j+1) - p(j-1))/min(p(j+1),p(j-1)) > parameter1,
	             u(j-1) - u(j+1)                     > 0,
		     e(j)/E(j)                           > parameter2 */
 
	for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	  for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	    for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
	      index = i + j*GridDimension[0] +
		      k*GridDimension[1]*GridDimension[0];
 
	      DelPressure = fabs(
		      (Pressure[index + Offset] -
		       Pressure[index - Offset]  ) /
		   min(Pressure[index + Offset],
		       Pressure[index - Offset]  )  );
 
	      DelVelocity = BaryonField[Vel1Num+dim][index - Offset] -
		            BaryonField[Vel1Num+dim][index + Offset];
 
	      EnergyThisCell = (Pressure[index])/(Gamma-1.0);
 
	      Energy1        = BaryonField[TENum  ][index + Offset] *
		               BaryonField[DensNum][index + Offset] ;
	      Energy2        = BaryonField[TENum  ][index         ] *
			       BaryonField[DensNum][index         ] ;
	      Energy3        = BaryonField[TENum  ][index - Offset] *
			       BaryonField[DensNum][index - Offset] ;
 
	      EnergyRatio = EnergyThisCell/max(max(Energy1, Energy2), Energy3);
 
//	      EnergyRatio =  Pressure[index]/
//		            (BaryonField[DensNum][index]*(Gamma-1.0)*
//			     BaryonField[TENum  ][index]);
 
	      FlaggingField[index] +=
		(DelPressure > MinimumPressureJumpForRefinement &&
		 DelVelocity > 0                                &&
		 EnergyRatio > MinimumEnergyRatioForRefinement     ) ? 1 : 0;
	    }
 
	Offset *= GridDimension[dim];
 
      }  // end loop over dimension
 
  /* clean up */
 
  delete [] Pressure;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
 
}
