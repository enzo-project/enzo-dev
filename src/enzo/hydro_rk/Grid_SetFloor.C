/***********************************************************************
/
/  GRID CLASS (SET INTERNAL ENERGY FLOOR)
/
/  written by: Peng Wang
/  date:       October, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::SetFloor()
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, 1.0);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, HMNum, H2INum, H2IINum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

#if 0
  float vx, vy, vz, v2, eint, emin, rho;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

        int igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
        rho = BaryonField[iden][igrid];
        vx = BaryonField[ivx][igrid];
        vy = BaryonField[ivy][igrid];
        vz = BaryonField[ivz][igrid];
        v2 = vx*vx + vy*vy + vz*vz;
        if (DualEnergyFormalism) {
          eint = BaryonField[ieint][igrid];
        }
        else {
          eint = BaryonField[ietot][igrid] - 0.5*v2;
	}

        emin = 4.0*0.48999*rho*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));

        eint = max(eint, emin);

        BaryonField[ietot][igrid] = eint + 0.5*v2;

        if (DualEnergyFormalism) {
          BaryonField[ieint][igrid] = eint;
        }

      } 
    }
  }
#endif

  if (HydroMethod == MHD_RK) {
    const float ca_min = MaximumAlvenSpeed;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  
	  int igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
	  float rho = BaryonField[DensNum][igrid];
	  float Bx = BaryonField[B1Num][igrid];
	  float By = BaryonField[B2Num][igrid];
	  float Bz = BaryonField[B3Num][igrid];
	  
	  float B2 = Bx*Bx + By*By + Bz*Bz;
	  float ca = sqrt(B2)/sqrt(rho);
	  
	  if (ca > ca_min) {
	    BaryonField[TENum][igrid] -= 0.5*B2/rho;
	    float rho1 = B2/pow(ca_min,2);
	    BaryonField[DensNum][igrid] = rho1;
	    BaryonField[TENum][igrid] += 0.5*B2/rho1;
	    printf("floor set based on MaximumAlvenSpeed: (%"GSYM" %"GSYM" %"GSYM"), rho: %"GSYM"->%"GSYM"\n", CellLeftEdge[0][i],
		   CellLeftEdge[1][j], CellLeftEdge[2][k], rho*DensityUnits, rho1*DensityUnits);
	  }
	}
      }
    }
  }
  
  return SUCCESS;
}
