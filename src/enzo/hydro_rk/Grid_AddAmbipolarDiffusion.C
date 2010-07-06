/***********************************************************************
/
/  GRID CLASS (ADD AMBIPOLAR DIFFUSION TERM)
/
/  written by: Peng Wang
/  date:       July, 2008
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
#include "fortran.def"
#include "CosmologyParameters.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddAmbipolarDiffusion()
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
				   TENum, B1Num, B2Num, B3Num, PhiNum);

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  double A = 1.4/(3.5e13*3e-16);
  A = A/(TimeUnits*sqrt(DensityUnits));

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* D = A/rho^1.5 * curlB x B x B. */

  double *D[3];
  for (int dim = 0; dim < 3; dim++) {
    D[dim] = new double[size];
  }

  /* Compute AD time step */
 
  FLOAT dt_ad = dtFixed, dt1;
  FLOAT dx = CellWidth[0][0];

  double rho, Bx, By, Bz, B2;
  int nx = (GridRank > 1) ? GridDimension[0] : 0, ny = (GridRank > 2) ? GridDimension[1] : 0;
  int igrid;
  int kmp = (GridRank > 2) ? 1 : 0,
    jmp = (GridRank > 1) ? 1 : 0;
  for (int k = GridStartIndex[2] - kmp; k <= GridEndIndex[2] + kmp; k++) {
    for (int j = GridStartIndex[1] - jmp; j <= GridEndIndex[1] + jmp; j++) {
      for (int i = GridStartIndex[0]-1; i <= GridEndIndex[0]+1; i++) {
	
	igrid = i + (j + k*ny)*nx;

	rho = BaryonField[DensNum][igrid];
	Bx  = BaryonField[B1Num][igrid];
	By  = BaryonField[B2Num][igrid];
	Bz  = BaryonField[B3Num][igrid];
	B2  = Bx*Bx + By*By + Bz*Bz;

	dt1 = 0.3*dx*pow(rho,1.5)/(A*B2)*dx;

	if (dt1 < dt_ad) dt_ad = dt1;

      }
    }
  }

  printf("AD time step = %"GSYM"\n", dt_ad);
  
  FLOAT dt_total = 0;
  while (dt_total < dtFixed) {
  
    /* First calculate D */
    
    double dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy;
    for (int k = GridStartIndex[2]-kmp; k <= GridEndIndex[2]+kmp; k++) {
      for (int j = GridStartIndex[1]-jmp; j <= GridEndIndex[1]+jmp; j++) {
	for (int i = GridStartIndex[0]-1; i <= GridEndIndex[0]+1; i++) {

	  igrid = i + (j + k*ny)*nx;
	  
	  rho = BaryonField[DensNum][igrid];
	  Bx  = BaryonField[B1Num][igrid];
	  By  = BaryonField[B2Num][igrid];
	  Bz  = BaryonField[B3Num][igrid];
	

	  dBxdy = (BaryonField[B1Num][igrid+nx   ] - BaryonField[B1Num][igrid-nx   ]);
	  dBxdz = (BaryonField[B1Num][igrid+nx*ny] - BaryonField[B1Num][igrid-nx*ny]);
	  dBydx = (BaryonField[B2Num][igrid+1    ] - BaryonField[B2Num][igrid-1    ]);
	  dBydz = (BaryonField[B2Num][igrid+nx*ny] - BaryonField[B2Num][igrid-nx*ny]);
	  dBzdx = (BaryonField[B3Num][igrid+1    ] - BaryonField[B3Num][igrid-1    ]);
	  dBzdy = (BaryonField[B3Num][igrid+nx   ] - BaryonField[B3Num][igrid-nx   ]);
	  
	  D[0][igrid] = Bx*Bz*(dBydx-dBxdy)-(By*By+Bz*Bz)*(dBzdy-dBydz)+Bx*By*(dBxdz-dBzdx);
	  D[1][igrid] = Bx*By*(dBzdy-dBydz)-(Bx*Bx+Bz*Bz)*(dBxdz-dBzdx)+By*Bz*(dBydx-dBxdy);
	  D[2][igrid] = By*Bz*(dBxdz-dBzdx)-(Bx*Bx+By*By)*(dBydx-dBxdy)+Bx*Bz*(dBzdy-dBydz);

	  D[0][igrid] *= A/dx/pow(rho,1.5);
	  D[1][igrid] *= A/dx/pow(rho,1.5);
	  D[2][igrid] *= A/dx/pow(rho,1.5);

       
	}
      }
    }
    
    /* Then calculte the AD term and update B field */

    double AD[3];
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  
	  igrid = i + (j + k*ny)*nx;
	  
	  AD[0] = (D[2][igrid+nx   ]-D[2][igrid-nx   ])-(D[1][igrid+nx*ny]-D[1][igrid-nx*ny]);
	  AD[1] = (D[0][igrid+nx*ny]-D[0][igrid-nx*ny])-(D[2][igrid+1    ]-D[2][igrid-1    ]);
	  AD[2] = (D[1][igrid+1    ]-D[1][igrid-1    ])-(D[0][igrid+nx   ]-D[0][igrid-nx   ]);

	  BaryonField[B1Num][igrid] += dt_ad/dx*AD[0];
	  BaryonField[B2Num][igrid] += dt_ad/dx*AD[1]; 
	  BaryonField[B3Num][igrid] += dt_ad/dx*AD[2];

	}
      }
    }
    
    dt_total += dt_ad;
    if (dt_total + dt_ad > dtFixed) dt_ad = dtFixed - dt_total;
    
  }

  /*for (int i = 0; i < size; i++) {
    printf("%"GSYM" ", BaryonField[B2Num][i]);
  }
  printf("\n");*/


  for (int dim = 0; dim < 3; dim++) {
    delete D[dim];
  }
  
  return SUCCESS;

}
