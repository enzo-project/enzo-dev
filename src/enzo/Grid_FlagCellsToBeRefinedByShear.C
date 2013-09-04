/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE PRESENCE OF SHEAR)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Alexei Kritsuk, Aug.2004
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
         float *TemperatureUnits, float *TimeUnits,
         float *VelocityUnits, FLOAT Time);
 
int grid::FlagCellsToBeRefinedByShear()
{
  /* declarations */
 
  int i, j, k, index, dim;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate temp array and fill with zeroes */
 
  float *DelVelocity = new float[size];
  for (i = 0; i < size; i++)
    DelVelocity[i] = 0.0;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  float c_sound2;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

/* Create temperature fields for calculating local sound speed */

  float *temperature = new float[size];

  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  /* loop over active dimensions */
 
  int Offset = 1;
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* For shear: [(du/dy)^2 + (dv/dx)^2 + (dw/dy)^2 +
                     (du/dz)^2 + (dv/dz)^2 + (dw/dx)^2  ] > parameter,
	 where:  du/dy = [u(j-1) - u(j+1)]/[2dy],
	 assume: dx=dy=dz=1;
	 parameter = epsilon_threshold * (c_sound/dx)^2
     */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
	    index = i + j*GridDimension[0] +
	                k*GridDimension[1]*GridDimension[0];
 
	    float DelVel1 = 0.0, DelVel2 = 0.0;
	    switch (dim) {
	    case 0:
	      if (GridRank > 1)
		DelVel1 = BaryonField[Vel1Num+1][index - Offset] -
		          BaryonField[Vel1Num+1][index + Offset];
	      if (GridRank > 2)
		DelVel2 = BaryonField[Vel1Num+2][index - Offset] -
		          BaryonField[Vel1Num+2][index + Offset];
	      break;
 
	    case 1:
	      DelVel1 = BaryonField[Vel1Num][index - Offset] -
	   	            BaryonField[Vel1Num][index + Offset];
	      if (GridRank > 2)
		  DelVel2 = BaryonField[Vel1Num+2][index - Offset] -
		            BaryonField[Vel1Num+2][index + Offset];
	      break;
 
	    case 2:
	      if (GridRank > 1)
		  DelVel1 = BaryonField[Vel1Num+1][index - Offset] -
		            BaryonField[Vel1Num+1][index + Offset];
	      DelVel2 = BaryonField[Vel1Num][index - Offset] -
		            BaryonField[Vel1Num][index + Offset];
	      break;
 
	    default:
	      break;
	    }
 
	    DelVel1 *= DelVel1;
	    DelVel2 *= DelVel2;
	    DelVelocity[index] += DelVel1 + DelVel2;
        // If ideal gas, then calculate local sound speed
        if (EOSType == 0) {
          c_sound2 = Gamma*kboltz*temperature[index]/(Mu*mh) / 
                     (VelocityUnits*VelocityUnits);
          //fprintf(stderr,"DelVelocity = %f\n", DelVelocity[index]);
          //fprintf(stderr,"epsilon*c_sound2 = %f\n", MinimumShearForRefinement*c_sound2);
        }
        // If nonideal gas, then assume sound speed is 1, as was
        // previously done for all calculations
        else {
          c_sound2 = 1.0;
        }
	    if (dim == GridRank-1) {
	      FlaggingField[index] +=
		  (DelVelocity[index] > (MinimumShearForRefinement*c_sound2)) ? 1 : 0;
        }
	  }
 
      Offset *= GridDimension[dim];
 
    }  // end loop over dimension
 
  /* clean up */
 
  delete DelVelocity;
  delete temperature;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i])

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
 
}
