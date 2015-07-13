/***********************************************************************
/
/  GRID CLASS (ADD VISCOSITY TERM)
/
/  written by: Peng Wang
/  date:       July, 2008
/  modified1:  11/09  Tom Abel, added ViscosityCoefficient
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

int grid::AddViscosity()
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

  double *d2Vx[5], *d2Vy[5], *d2Vz[5];
  for (int i = 0; i < 5; i++) {
    d2Vx[i] = new double[activesize];
    d2Vy[i] = new double[activesize];
    d2Vz[i] = new double[activesize];
  }

  /* Compute the kinematic viscosity */

  float *viscosity = new float[activesize];

  this->ComputeViscosity(viscosity, DensNum);

  /* Compute viscous time step */
 
  FLOAT dt_vis = dtFixed, dt1;
  FLOAT dx = CellWidth[0][0];

  for (int i = 0; i < activesize; i++) {
    dt1 = 0.1*dx/viscosity[i]*dx;
    if (dt1 < dt_vis) dt_vis = dt1;
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
	  d2Vx[0][n] = (BaryonField[Vel1Num][ip1] - 2.0*BaryonField[Vel1Num][igrid] + BaryonField[Vel1Num][im1]);
	  d2Vy[0][n] = (BaryonField[Vel2Num][ip1] - 2.0*BaryonField[Vel2Num][igrid] + BaryonField[Vel2Num][im1]);
	  d2Vz[0][n] = (BaryonField[Vel3Num][ip1] - 2.0*BaryonField[Vel3Num][igrid] + BaryonField[Vel3Num][im1]);
	  
	  int jp1 = i + (j+1+k*GridDimension[1])*GridDimension[0];
	  int jm1 = i + (j-1+k*GridDimension[1])*GridDimension[0];
	  d2Vx[1][n] = (BaryonField[Vel1Num][jp1] - 2.0*BaryonField[Vel1Num][igrid] + BaryonField[Vel1Num][jm1]);
	  d2Vy[1][n] = (BaryonField[Vel2Num][jp1] - 2.0*BaryonField[Vel2Num][igrid] + BaryonField[Vel2Num][jm1]);
	  d2Vz[1][n] = (BaryonField[Vel3Num][jp1] - 2.0*BaryonField[Vel3Num][igrid] + BaryonField[Vel3Num][jm1]);
	  
	  int kp1 = i + (j+(k+1)*GridDimension[1])*GridDimension[0];
	  int km1 = i + (j+(k-1)*GridDimension[1])*GridDimension[0];
	  d2Vx[2][n] = (BaryonField[Vel1Num][kp1] - 2.0*BaryonField[Vel1Num][igrid] + BaryonField[Vel1Num][km1]);
	  d2Vy[2][n] = (BaryonField[Vel2Num][kp1] - 2.0*BaryonField[Vel2Num][igrid] + BaryonField[Vel2Num][km1]);
	  d2Vz[2][n] = (BaryonField[Vel3Num][kp1] - 2.0*BaryonField[Vel3Num][igrid] + BaryonField[Vel3Num][km1]);
	  
	  int ip1jp1 = i+1 + (j+1+k*GridDimension[1])*GridDimension[0];
	  int ip1jm1 = i+1 + (j-1+k*GridDimension[1])*GridDimension[0];
	  int im1jp1 = i-1 + (j+1+k*GridDimension[1])*GridDimension[0];
	  int im1jm1 = i-1 + (j-1+k*GridDimension[1])*GridDimension[0];
	  d2Vx[3][n] = (BaryonField[Vel1Num][ip1jp1] + BaryonField[Vel1Num][im1jm1] -
			BaryonField[Vel1Num][ip1jm1] - BaryonField[Vel1Num][im1jp1]);
	  d2Vy[3][n] = (BaryonField[Vel2Num][ip1jp1] + BaryonField[Vel2Num][im1jm1] -
			BaryonField[Vel2Num][ip1jm1] - BaryonField[Vel2Num][im1jp1]);
	  
	  int ip1kp1 = i+1 + (j+(k+1)*GridDimension[1])*GridDimension[0];
	  int ip1km1 = i+1 + (j+(k-1)*GridDimension[1])*GridDimension[0];
	  int im1kp1 = i-1 + (j+(k+1)*GridDimension[1])*GridDimension[0];
	  int im1km1 = i-1 + (j+(k-1)*GridDimension[1])*GridDimension[0];
	  d2Vx[4][n] = (BaryonField[Vel1Num][ip1kp1] + BaryonField[Vel1Num][im1km1] -
			BaryonField[Vel1Num][ip1km1] - BaryonField[Vel1Num][im1kp1]);
	  d2Vz[3][n] = (BaryonField[Vel3Num][ip1kp1] + BaryonField[Vel3Num][im1km1] -
			BaryonField[Vel3Num][ip1km1] - BaryonField[Vel3Num][im1kp1]);
	  
	  int jp1kp1 = i + (j+1+(k+1)*GridDimension[1])*GridDimension[0];
	  int jp1km1 = i + (j+1+(k-1)*GridDimension[1])*GridDimension[0];
	  int jm1kp1 = i + (j-1+(k+1)*GridDimension[1])*GridDimension[0];
	  int jm1km1 = i + (j-1+(k-1)*GridDimension[1])*GridDimension[0];
	  d2Vy[4][n] = (BaryonField[Vel2Num][jp1kp1] + BaryonField[Vel2Num][jm1km1] -
			BaryonField[Vel2Num][jp1km1] - BaryonField[Vel2Num][jm1kp1]);
	  d2Vz[4][n] = (BaryonField[Vel3Num][jp1kp1] + BaryonField[Vel3Num][jm1km1] -
			BaryonField[Vel3Num][jp1km1] - BaryonField[Vel3Num][jm1kp1]);
	  
	}
      }
    }
    
    /* Update velocity */
    
    n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i + (j+k*GridDimension[1])*GridDimension[0];
	  
	  BaryonField[Vel1Num][igrid] += 
	    dt_vis/dx*viscosity[n]/dx*
	    ((d2Vx[0][n]+d2Vx[1][n]+d2Vx[2][n]) + 1.0/3.0*(d2Vx[0][n]+d2Vy[3][n]+d2Vz[3][n]));
	  
	  BaryonField[Vel2Num][igrid] += 
	    dt_vis/dx*viscosity[n]/dx*
	    ((d2Vy[0][n]+d2Vy[1][n]+d2Vy[2][n]) + 1.0/3.0*(d2Vx[3][n]+d2Vy[1][n]+d2Vz[4][n]));
	  
	  BaryonField[Vel3Num][igrid] +=
	    dt_vis/dx*viscosity[n]/dx*
	    ((d2Vz[0][n]+d2Vz[1][n]+d2Vz[2][n]) + 1.0/3.0*(d2Vx[4][n]+d2Vy[4][n]+d2Vz[2][n]));
	  
	}
      }
    }
    
    dt_total += dt_vis;
    if (dt_total + dt_vis > dtFixed) dt_vis = dtFixed - dt_total;
  }

  for (int i = 0; i < 5; i++) {
    delete d2Vx[i];
    delete d2Vy[i];
    delete d2Vz[i];
  }
  
  delete viscosity;
  
  return SUCCESS;

}

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::ComputeViscosity(float *viscosity, int DensNum)
{
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  int n = 0, igrid;  
  if (UseViscosity == 1) {
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) 
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) 
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) 
	  viscosity[n] = ViscosityCoefficient;
    //printf("VISC: %"FSYM"\n", ViscosityCoefficient);
  }
  else if (UseViscosity == 2)   {
    double alpha = 0.1;
    double Msun = 1.989e33;
    double GravConst = 6.672e-8;
    double M;
    if (ExternalGravity == 4) M = 2.0*ExternalGravityDensity*Msun;
    if (ExternalGravity == 5) M = ExternalGravityDensity*Msun;
    
    float Omega, rho, p, eint, h, cs, dpdrho, dpde;
    FLOAT x, y, R;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  
	  igrid = i + (j+k*GridDimension[1])*GridDimension[0];
	  
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];	
	  R = sqrt(pow(x-0.5,2) + pow(y-0.5,2));
	  R = max(R, 0.5*CellWidth[0][0]);
	  
	  rho = BaryonField[DensNum][igrid];
	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1);
	  
	  Omega = sqrt(GravConst*M/pow(R*LengthUnits,3))*TimeUnits;
	  
	  viscosity[n] = alpha*cs*cs/Omega;
	  
	}
      }
    }
  }
  
  return SUCCESS;
  
}
