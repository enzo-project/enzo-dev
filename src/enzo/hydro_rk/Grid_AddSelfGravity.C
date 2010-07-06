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
#include "TopGridData.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddSelfGravity(float coef)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  if (Coordinate == Cartesian) {
    int i, j, k, igrid;
    float temp1, temp2, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, igrid++) {
	  temp1 = OldBaryonField[DensNum][igrid] / BaryonField[DensNum][igrid];
	  temp2 =  dtFixed*0.5*(1.0 + temp1);
	  gx = AccelerationField[0][igrid];
	  gy = (GridRank > 1) ? AccelerationField[1][igrid] : 0.0;
	  gz = (GridRank > 2) ? AccelerationField[2][igrid] : 0.0;
	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];
	  vx_old = OldBaryonField[Vel1Num][igrid];
	  vy_old = OldBaryonField[Vel2Num][igrid];
	  vz_old = OldBaryonField[Vel3Num][igrid];

	  /*BaryonField[Vel1Num][igrid] += coef*gx*temp1*dtFixed;
	  BaryonField[Vel2Num][igrid] += coef*gy*temp1*dtFixed;
	  BaryonField[Vel3Num][igrid] += coef*gz*temp1*dtFixed;
	  BaryonField[TENum][igrid] += ceof**/

	  BaryonField[Vel1Num][igrid] += coef*gx*temp2;
	  BaryonField[Vel2Num][igrid] += coef*gy*temp2;
	  BaryonField[Vel3Num][igrid] += coef*gz*temp2;

	  BaryonField[TENum][igrid] += 
	    coef*dtFixed*0.5*(gx*(vx+vx_old*temp1)+gy*(vy+vy_old*temp1)+gz*(vz+vz_old*temp1));
	}
      }
    }
  }


  return SUCCESS;
}
