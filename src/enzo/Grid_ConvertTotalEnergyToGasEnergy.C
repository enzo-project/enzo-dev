/***********************************************************************
/
/  GRID CLASS (CONVERT TOTAL ENERGY TO GAS ENERGY)
/
/  written by: Greg Bryan
/  date:       February, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int grid::ConvertTotalEnergyToGasEnergy()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  int dim, i, size = 1;
 
  /* Compute the size of the grid. */
 
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "CTETGE: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* Subtract kinetic component. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[TENum][i] -=
	0.5*BaryonField[Vel1Num+dim][i]*BaryonField[Vel1Num+dim][i];
 
  return SUCCESS;
}
