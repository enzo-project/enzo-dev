/***********************************************************************
/
/  GRID CLASS (ADD ACCELERATION FROM RADIATION PRESSURE)
/
/  written by: John H. Wise
/  date:       December, 2005
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
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

int grid::AddRadiationPressureAcceleration()
{

  /* Return if this does not concern us */
  if (!(RadiativeTransfer && RadiationPressure)) return SUCCESS;

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int dim, i, j, k, index, size = 1;

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

  /* Get acceleration fields from radiation pressure */

  int RPresNum1, RPresNum2, RPresNum3;
  if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) 
      == FAIL) {
    ENZO_FAIL("Error in IdentifyRadiationPressureFields.\n");
  }

  /* Add acceleration fields from radiation pressure */
  index = 0;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	for (dim = 0; dim < GridRank; dim++) {
	  AccelerationField[dim][index] += BaryonField[RPresNum1+dim][index];
	  /*
	  if (fabs(BaryonField[RPresNum1+dim][index]) > 

	      fabs(0.05*AccelerationField[dim][index]))  
	    fprintf(stdout, "AddRPAccel[dim %"ISYM" :: %"ISYM" %"ISYM" %"ISYM"]: "
		    "Accel = %"GSYM", RPAccel = %"GSYM"\n", 
		    dim, i, j, k, AccelerationField[dim][index],
		    BaryonField[RPresNum1+dim][index]);  
	  */
	}
    }  // ENDFOR j


  return SUCCESS;
}

