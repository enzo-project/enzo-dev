/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE PRESENCE OF SHEAR)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Alexei Kritsuk, Aug.2004
/  modified2:  Cameron Hummels, Oct.2013
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
 
int grid::FlagCellsToBeRefinedByShear(int level)
{
  /* declarations */
 
  int i, j, k, index, dim, LeftOffset, RightOffset;
  float Shear2, Denom, Norm2;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate temp arrays and fill with zeroes */
 
  float *DelVel1= new float[size];
  float *DelVel2 = new float[size];
  float *DelVel3 = new float[size];

  for (i = 0; i < size; i++) {
    DelVel1[i] = 0.0;
    DelVel2[i] = 0.0;
    DelVel3[i] = 0.0;
  }

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  float c_sound2;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

/* Create temperature fields for calculating local sound speed */

  float *pressure = new float[size];

  if (this->ComputePressure(Time, pressure) == FAIL){
    ENZO_FAIL("Error in grid->ComputePressure.");
  }

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  /* Set up stencil for Zeus vs PPM data, since PPM uses cell-centered 
     velocities and Zeus uses edge-centered velocities */
  if (HydroMethod == 2) { // Zeus
    LeftOffset = 1;
    RightOffset = 0;
    Denom = 1.0;
  }
  else {                 // PPM
    LeftOffset = 1;
    RightOffset = 1;
    Denom = 2.0;
  }

  /* loop over active dimensions */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* Absolute Shear: [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 + 
                          (dvx/dz + dvz/dx)^2  ]^(0.5)
	     where:  dvx/dy = [vx(j-1) - vx(j+1)]/[2dy],
         in units of s^(-1)
    
         We scale this to be dimensionless by multiplying by (dx/c_sound),
         which is equivalent to dividing out all of the dx/dy/dz values from 
         the shear equation and dividing by c_sound.  This yields a 
         discontinuous quantity, so we scale it by 2**level, which means
         the shear is relative to the root grid dx/dy/dz and the local sound
         speed:

      Dimensionless Shear: [(dvx + dvy)^2 + (dvz + dvy)^2 + 
                            (dvx + dvz)^2  ]^(0.5) * (c_sound)
    
       
      Refinement occurs if: Dimensionless Shear > epsilon_threshold
      where epsilon_threshold is a user-defined value tied to the local mach
      number.

      N.B.: DelVel1 = dvx/dy + dvy/dx
            DelVel2 = dvx/dz + dvz/dx
            DelVel3 = dvz/dy + dvy/dz
     */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	  for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
	    index = i + j*GridDimension[0] +
	                k*GridDimension[1]*GridDimension[0];
 
	    switch (dim) {
	    case 0: 
	      if (GridRank > 1) // dvy/dx
          DelVel1[index] += (BaryonField[Vel2Num][index - LeftOffset] -
                             BaryonField[Vel2Num][index + RightOffset]) /
                             this->CellWidth[dim][k];
	      if (GridRank > 2) // dvz/dx
          DelVel2[index] += (BaryonField[Vel3Num][index - LeftOffset] -
                             BaryonField[Vel3Num][index + RightOffset]) /
                             this->CellWidth[dim][k];
	      break;
 
	    case 1:             // dvx/dy
	      DelVel1[index] += (BaryonField[Vel1Num][index - LeftOffset] -
	   	                     BaryonField[Vel1Num][index + RightOffset]) /
                             this->CellWidth[dim][k];
	      if (GridRank > 2) // dvz/dy
		  DelVel3[index] += (BaryonField[Vel3Num][index - LeftOffset] -
		                     BaryonField[Vel3Num][index + RightOffset]) /
                             this->CellWidth[dim][k];
	      break;
 
	    case 2:             // dvx/dz and dvy/dz 
	      DelVel2[index] += (BaryonField[Vel1Num][index - LeftOffset] -
		                     BaryonField[Vel1Num][index + RightOffset]) /
                             this->CellWidth[dim][k];
		  DelVel3[index] += (BaryonField[Vel2Num][index - LeftOffset] -
		                     BaryonField[Vel2Num][index + RightOffset]) / 
                             this->CellWidth[dim][k];
	      break;
 
	    default:
	      break;
	    }
 
	    if (dim == GridRank-1) {    // if we're done

        // If ideal gas, then calculate local sound speed
        if (EOSType == 0) 
          c_sound2 = Gamma*pressure[index]/BaryonField[DensNum][index];
        
        // If nonideal gas, then assume sound speed is 1
        else 
          c_sound2 = 1.0;
        
        // Calculate the normalization to make sure shear is continuous
        // across refinement boundaries
        // Norm2 = pow(2.0, level) * pow(2.0, level);

        // check to make sure width decreases by a factor of 2 when level increases by 1
        // fprintf(stderr,"Level = %i; k Width = %g; j Width = %g; i width = %g\n", 
        //         level, this->CellWidth[dim][k], this->CellWidth[dim][j], this->CellWidth[dim][i]);

        Shear2 = (DelVel1[index]*DelVel1[index] + 
                  DelVel2[index]*DelVel2[index] + 
                  DelVel3[index]*DelVel3[index]) / 
                 (Denom*Denom*c_sound2);

        FlaggingField[index] +=
        (Shear2 > (MinimumShearForRefinement*MinimumShearForRefinement) / 
                  (this->CellWidth[dim][k]*this->CellWidth[dim][k]) ) ? 1 : 0;
        }
	  }
 
      LeftOffset *= GridDimension[dim];
      RightOffset *= GridDimension[dim];
 
    }  // end loop over dimension
 
  /* clean up */
 
  delete [] DelVel1;
  delete [] DelVel2;
  delete [] DelVel3;
  delete [] pressure;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i])

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
 
}
