/***********************************************************************
/
/  GRID CLASS (RESET MAGNETIC FIELD)
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1:  

/  PURPOSE: 
/
/  RETURNS:
/    SUCCESS or FAIL
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
#include "fortran.def"
#include "CosmologyParameters.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
  
int grid::MagneticFieldResetter(int level)
{

  if (ResetMagneticField == FALSE)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = NumberOfGhostZones;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;

  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
 
  /* Compute the redshift. */
 
  float zred;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  zred = 1.0*(1.0+InitialRedshift)/a - 1.0;
 
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;

  /* Convert magnetic amplitude into code unit */

  for (dim = 0; dim < GridRank; dim++)
    ResetMagneticFieldAmplitude[dim] /= MagneticUnits;

  /* Reset the magnetic field values, update the total energy */

  for (int i = 0; i < size; i++) {

    /* Set the magnetic field value when the density is bigger than the critical density */

    if (BaryonField[DensNum][i] > 1.0) {   

      if(DualEnergyFormalism)
	BaryonField[TENum][i] -= 0.5*(BaryonField[B1Num][i]*BaryonField[B1Num][i] + 
				      BaryonField[B1Num][i]*BaryonField[B1Num][i] + 
				      BaryonField[B1Num][i]*BaryonField[B1Num][i])/BaryonField[DensNum][i];
      
      BaryonField[B1Num][i]  = ResetMagneticFieldAmplitude[0];
      BaryonField[B2Num][i]  = ResetMagneticFieldAmplitude[1];
      BaryonField[B3Num][i]  = ResetMagneticFieldAmplitude[2];
      BaryonField[PhiNum][i] = 0.0;
      
      if(DualEnergyFormalism)
	BaryonField[TENum][i] += 0.5*(BaryonField[B1Num][i]*BaryonField[B1Num][i] + 
				      BaryonField[B1Num][i]*BaryonField[B1Num][i] + 
				      BaryonField[B1Num][i]*BaryonField[B1Num][i])/BaryonField[DensNum][i];
    } else {

      BaryonField[B1Num][i]  = 0.0;
      BaryonField[B2Num][i]  = 0.0;
      BaryonField[B3Num][i]  = 0.0;
      BaryonField[PhiNum][i] = 0.0;

    }
      
  }

  /* For future use, convert the amplitude back into Gauss */

  for (dim = 0; dim < GridRank; dim++)
    ResetMagneticFieldAmplitude[dim] *= MagneticUnits;

  return SUCCESS;

}
