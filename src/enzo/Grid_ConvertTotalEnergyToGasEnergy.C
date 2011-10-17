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
#include "ErrorExceptions.h"
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
    ENZO_FAIL("CTETGE: Error in IdentifyPhysicalQuantities.\n");
  }
  printf("DensNum = %d, GENum = %d, vel1Num = %d, vel2Num = %d,vel3Num = %d,TENum = %d \n",DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum ); 
  /* Subtract kinetic component. */
  //printf("BaryonField[DensNum][1] = %g \n",i,BaryonField[DensNum][1]); 
  //printf("BaryonField[TENum][2] = %g \n",i,BaryonField[TENum][2]); 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++){
      //printf("BaryonField[TENum][2] = %g, BaryonField[Vel1Num] = %g, BaryonField[Vel2Num] = %g, BaryonField[Vel3Num] = %g  \n",BaryonField[TENum][2], BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num]); 

      BaryonField[TENum][i] -= 0.5*BaryonField[Vel1Num+dim][i]*BaryonField[Vel1Num+dim][i];
      if(BaryonField[TENum][i] < 0.0) 

	printf(" PROBLEM!!!! BaryonField[TENum][i] = %g, BaryonField[Vel1Num] = %g, BaryonField[Vel2Num] = %g, BaryonField[Vel3Num] = %g  \n",BaryonField[TENum][i], BaryonField[Vel1Num][i], BaryonField[Vel2Num][i], BaryonField[Vel3Num][i]);
    }
  return SUCCESS;
}
