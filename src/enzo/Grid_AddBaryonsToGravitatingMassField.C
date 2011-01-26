/***********************************************************************
/
/  GRID CLASS (ADD THE BARYON MASS TO THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
 
int FindField(int f, int farray[], int n);
 
#define INCLUDE_GHOST_ZONES
 
int grid::AddBaryonsToGravitatingMassField()
{
  /* Don't bother if there are no baryon fields or if this is not the
     grid's own processor. */
 
  if (NumberOfBaryonFields == 0 || MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  int dim, i, j, k;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
  }
 
  /* Error check. */
 
  if (GravitatingMassField == NULL) {
    ENZO_FAIL("GravitatingMassField not allocated.\n");
  }
 
  /* Compute Offset between baryon field and GravitatingMassField. */
 
  int Offset[MAX_DIMENSION] = {0, 0, 0};
  for (dim = 0; dim < GridRank; dim++) {
    Offset[dim] = nint((CellLeftEdge[dim][0] -
			GravitatingMassFieldLeftEdge[dim])/CellWidth[dim][0]);
    if (Offset[dim] < 0) {
      ENZO_VFAIL("Offset[%"ISYM"] = %"ISYM" < 0\n", dim, Offset[dim])

    }
  }
 
  /* Loop over mesh, adding to field. */
 
  int gmfindex, index = 0;
 
#ifdef INCLUDE_GHOST_ZONES
 
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++) {
      gmfindex = ((k+Offset[2])*GravitatingMassFieldDimension[1] +
		  (j+Offset[1]))*GravitatingMassFieldDimension[0] +
		     Offset[0];
      for (i = 0; i < GridDimension[0]; i++, index++, gmfindex++)
	GravitatingMassField[gmfindex] += BaryonField[DensNum][index];
    }
 
#else  /* INCLUDE_GHOST_ZONES */
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      gmfindex = ((k+Offset[2])*GravitatingMassFieldDimension[1] +
		  (j+Offset[1]))*GravitatingMassFieldDimension[0] +
		   GridStartIndex[0]+Offset[0];
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++,
	     index++, gmfindex++)
	GravitatingMassField[gmfindex] += BaryonField[DensNum][index];
    }
 
#endif /* INCLUDE_GHOST_ZONES */
 
  return SUCCESS;
}
