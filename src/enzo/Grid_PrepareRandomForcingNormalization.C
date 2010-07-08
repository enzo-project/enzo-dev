/***********************************************************************
/
/  GRID CLASS (PREPARE NORMALIZATION FOR RANDOM FORCING FIELDS +
/              COMPUTE GLOBAL VALUES OF INTEREST)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::PrepareRandomForcingNormalization(float * GlobVal, int GlobNum)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("GPRFN: Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Allocate field and compute temperature (it is actually c^2) */
 
  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  float *temperature = new float[size];
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");

  }
 
  /* Loop over active part of fields and sum up each of the required
     quantities: rho*v*(delta v) and rho*(delta v)^2; we also want
      to compute here rho*v^2/c^2 and v^2/c^2 as we will need M_rms
     (mass and volume weighted) for the level 0 grid(s) (only!!).
     Also note: temperature is actually c^2. */
 
  for (int num = 0; num < GlobNum; num++)
    GlobVal[num] = 0.0;
 
  GlobVal[7] = huge_number; // will be used for min density
 
  int i, j, k, index;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	for (dim = 0; dim < GridRank; dim++) {
	  int vel = Vel1Num + dim;
	  GlobVal[0] += BaryonField[vel][index]*
		        RandomForcingField[dim][index]*
	                BaryonField[DensNum][index];
	  GlobVal[1] += RandomForcingField[dim][index]*
		        RandomForcingField[dim][index]*
	                BaryonField[DensNum][index];
	  GlobVal[2] += BaryonField[vel][index]*
		        BaryonField[vel][index]*
	                BaryonField[DensNum][index]/
	                temperature[index];
	  GlobVal[3] += BaryonField[vel][index]*
		        BaryonField[vel][index]/
	                temperature[index];
	  GlobVal[4] += BaryonField[vel][index]*
		        BaryonField[vel][index]*
	                BaryonField[DensNum][index];
	  GlobVal[5] += BaryonField[vel][index]*
		        BaryonField[vel][index];
	}
	GlobVal[6] += BaryonField[DensNum][index]*
	  BaryonField[DensNum][index];
	GlobVal[7]  = min(GlobVal[7], BaryonField[DensNum][index]);
	GlobVal[8]  = max(GlobVal[8], BaryonField[DensNum][index]);
      }
 
  /* clean up. */
 
  delete temperature;
 
  return SUCCESS;
 
}
