/***********************************************************************
/
/  GRID CLASS (ADD RESISTIVITY TERM)
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


int grid::AddResistivity()
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


  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
  }

  double *d2Bx[3], *d2By[3], *d2Bz[3];
  for (int dim = 0; dim < 3; dim++) {
    d2Bx[dim] = new double[activesize];
    d2By[dim] = new double[activesize];
    d2Bz[dim] = new double[activesize];

    for (int i = 0; i < activesize; i++) {
      d2Bx[dim][i] = 0.0;
      d2By[dim][i] = 0.0;
      d2Bz[dim][i] = 0.0;
    }
  }

  /* Compute the kinematic viscosity */

  float *resistivity = new float[activesize];

  this->ComputeResistivity(resistivity, DensNum);

  /* Compute resistive time step */
 
  FLOAT dt_res = dtFixed, dt1;
  FLOAT dx = CellWidth[0][0];

  for (int i = 0; i < activesize; i++) {
    dt1 = 0.3*dx/resistivity[i]*dx;
    if (dt1 < dt_res) dt_res = dt1;
  }

  FLOAT dt_total = 0;

  while (dt_total < dtFixed) {

    /* Calculate the 2nd derivatives in the viscosity term */
  
    int n = 0, igrid;

    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  
	  igrid = i + (j+k*GridDimension[1])*GridDimension[0];
	  
	  int ip1 = i+1 + (j+k*GridDimension[1])*GridDimension[0];
	  int im1 = i-1 + (j+k*GridDimension[1])*GridDimension[0];
	  d2Bx[0][n] = (BaryonField[B1Num][ip1] - 2.0*BaryonField[B1Num][igrid] + BaryonField[B1Num][im1]);
	  d2By[0][n] = (BaryonField[B2Num][ip1] - 2.0*BaryonField[B2Num][igrid] + BaryonField[B2Num][im1]);
	  d2Bz[0][n] = (BaryonField[B3Num][ip1] - 2.0*BaryonField[B3Num][igrid] + BaryonField[B3Num][im1]);

	  if (GridRank > 1) {
	    int jp1 = i + (j+1+k*GridDimension[1])*GridDimension[0];
	    int jm1 = i + (j-1+k*GridDimension[1])*GridDimension[0];
	    d2Bx[1][n] = (BaryonField[B1Num][jp1] - 2.0*BaryonField[B1Num][igrid] + BaryonField[B1Num][jm1]);
	    d2By[1][n] = (BaryonField[B2Num][jp1] - 2.0*BaryonField[B2Num][igrid] + BaryonField[B2Num][jm1]);
	    d2Bz[1][n] = (BaryonField[B3Num][jp1] - 2.0*BaryonField[B3Num][igrid] + BaryonField[B3Num][jm1]);
	  }
	  
	  if (GridRank > 2) {
	    int kp1 = i + (j+(k+1)*GridDimension[1])*GridDimension[0];
	    int km1 = i + (j+(k-1)*GridDimension[1])*GridDimension[0];
	    d2Bx[2][n] = (BaryonField[B1Num][kp1] - 2.0*BaryonField[B1Num][igrid] + BaryonField[B1Num][km1]);
	    d2By[2][n] = (BaryonField[B2Num][kp1] - 2.0*BaryonField[B2Num][igrid] + BaryonField[B2Num][km1]);
	    d2Bz[2][n] = (BaryonField[B3Num][kp1] - 2.0*BaryonField[B3Num][igrid] + BaryonField[B3Num][km1]);
	  }
	  	  
	}
      }
    }
    
    /* Update Magnetic Field */
    
    n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i + (j+k*GridDimension[1])*GridDimension[0];

	  BaryonField[B1Num][igrid] += dt_res/dx*resistivity[n]/dx*(d2Bx[0][n]+d2Bx[1][n]+d2Bx[2][n]);	  
	  BaryonField[B2Num][igrid] += dt_res/dx*resistivity[n]/dx*(d2By[0][n]+d2By[1][n]+d2By[2][n]);	  
	  BaryonField[B3Num][igrid] += dt_res/dx*resistivity[n]/dx*(d2Bz[0][n]+d2Bz[1][n]+d2Bz[2][n]);
	  
	}
      }
    }
    
    dt_total += dt_res;
    if (dt_total + dt_res > dtFixed) dt_res = dtFixed - dt_total;

  }

  for (int dim = 0; dim < 3; dim++) {
    delete d2Bx[dim];
    delete d2By[dim];
    delete d2Bz[dim];
  }

  delete resistivity;
  
  return SUCCESS;

}

int grid::ComputeResistivity(float *resistivity, int DensNum)
{

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  double ResistivityUnits = pow(LengthUnits,2)/TimeUnits;

  int n = 0, igrid;
  float rho, p, eint, h, cs, dpdrho, dpde, T, nH2;
  double mp = 1.67e-24;
  FLOAT x, y, R;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	
	igrid = i + (j+k*GridDimension[1])*GridDimension[0];
	
	/* Calculate molecular number density */

	rho = BaryonField[DensNum][igrid];
	nH2 = rho*DensityUnits/(Mu*mp);

	/* Calculate temperature */

	EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1);
	T = p*Mu/rho*TemperatureUnits;
	
	resistivity[n] = 740.0/(5.7e-4)*nH2*sqrt(T/10.0)*(1.0-tanh(nH2/1e15));
	resistivity[n] /= ResistivityUnits;
	
      }
    }
  }

  return SUCCESS;

}
