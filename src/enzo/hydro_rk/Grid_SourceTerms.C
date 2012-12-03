/***********************************************************************
/
/  GRID CLASS (ADD SOURCE TERMS)
/
/  written by: Peng Wang
/  date:       May, 2007
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


int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::SourceTerms(float **dU)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, HMNum, H2INum, H2IINum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }


    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */

    FLOAT a = 1, dadt;
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	  == FAIL) {
	ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
      }


  if (DualEnergyFormalism) {   
    if (Coordinate == Cartesian) {
      int igrid, ip1, im1, jp1, jm1, kp1, km1;
      FLOAT dtdx = 0.5*dtFixed/CellWidth[0][0]/a,
	dtdy = (GridRank > 1) ? 0.5*dtFixed/CellWidth[1][0]/a : 0.0,
	dtdz = (GridRank > 2) ? 0.5*dtFixed/CellWidth[2][0]/a : 0.0;
      float rho, eint, p, divVdt, h, cs, dpdrho, dpde;
      float min_coeff = 0.0;
      if (UseMinimumPressureSupport) {
	min_coeff = MinimumPressureSupportParameter*
	  0.48999*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));
      }
      int n = 0;
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	    igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	    ip1 = igrid + 1;
	    im1 = igrid - 1;
	    jp1 = (GridRank > 1) ? i + GridDimension[0]*(j + 1 + k*GridDimension[1]) : 0;
	    jm1 = (GridRank > 1) ? i + GridDimension[0]*(j - 1 + k*GridDimension[1]) : 0;
	    kp1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k+1)*GridDimension[1]) : 0;
	    km1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k-1)*GridDimension[1]) : 0;
	    divVdt = dtdx*(BaryonField[Vel1Num][ip1] - BaryonField[Vel1Num][im1]) +
	      dtdy*(BaryonField[Vel2Num][jp1] - BaryonField[Vel2Num][jm1]) +
	      dtdz*(BaryonField[Vel3Num][kp1] - BaryonField[Vel3Num][km1]);
	    rho = BaryonField[DensNum][igrid];
	    eint = BaryonField[GENum][igrid];
	    eint = max(eint, min_coeff*rho);
	    EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
	    dU[iEint][n] -= p*divVdt;
	  }
	}
      }
    }

    if (Coordinate == Spherical) {
      int igrid, ip1, im1, jp1, jm1, kp1, km1;
      FLOAT xl, xc, xr, yl, yc, yr;
      FLOAT dx = CellWidth[0][0],
	dy = (GridRank > 1) ? CellWidth[1][0] : 1.0,
	dz = (GridRank > 1) ? CellWidth[2][0] : 1.0;
      FLOAT dtdx = 0.5*dtFixed/dx,
	dtdy = (GridRank > 1) ? 0.5*dtFixed/dy : 0.0,
	dtdz = (GridRank > 2) ? 0.5*dtFixed/dz : 0.0;
      float rho, eint, p, divVdt, h, cs, dpdrho, dpde;
      int n = 0;
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	    igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	    xl = CellLeftEdge[0][i] - 0.5*dx;
	    xc = xl + dx;
	    xr = xc + dx;
	    yl = (GridRank > 1) ? CellLeftEdge[1][j] - 0.5*dy : 0.0;
	    yc = yl + dy;
	    yr = yc + dy;
	    ip1 = igrid + 1;
	    im1 = igrid - 1;
	    jp1 = (GridRank > 1) ? i + GridDimension[0]*(j + 1 + k*GridDimension[1]) : 0;
	    jm1 = (GridRank > 1) ? i + GridDimension[0]*(j - 1 + k*GridDimension[1]) : 0;
	    kp1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k+1)*GridDimension[1]) : 0;
	    km1 = (GridRank > 2) ? i + GridDimension[0]*(j + (k-1)*GridDimension[1]) : 0;
	    divVdt = dtdx*(xr*xr*BaryonField[Vel1Num][ip1] - xl*xl*BaryonField[Vel1Num][im1])/(xc*xc) +
	      dtdy*(sin(yr)*BaryonField[Vel2Num][jp1] - sin(yl)*BaryonField[Vel2Num][jm1])/(xc*sin(yc)) +
	      dtdz*(BaryonField[Vel3Num][kp1] - BaryonField[Vel3Num][km1])/(xc*sin(yc));
	    rho = BaryonField[DensNum][igrid];
	    eint = BaryonField[GENum][igrid];
	    EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
	    dU[iEint][n] -= p*divVdt;
	  }
	}
      }
    }
  }

  if (Coordinate == Spherical) {
    float rho, etot, eint, vx, vy, vz, v2, e, h, cs, p, dpdrho, dpde, coty;
    FLOAT x, y, dtxinv;
    float pi = 4.0*atan(1.0);
    int n = 0, igrid;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = i+(j+k*GridDimension[1])*GridDimension[0];
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
          y = (GridRank > 1) ? CellLeftEdge[1][j] + 0.5*CellWidth[1][j] : pi*0.5;
          rho = BaryonField[DensNum ][igrid];
          vx  = BaryonField[Vel1Num  ][igrid];
          vy  = BaryonField[Vel2Num  ][igrid];
          vz  = BaryonField[Vel3Num  ][igrid];
	  if (DualEnergyFormalism) {
	    eint = BaryonField[GENum][igrid];
	  }
	  else {
	    etot = BaryonField[TENum][igrid];
	    v2 = vx*vx + vy*vy + vz*vz;
	    eint = etot - 0.5*v2;
	  }
                  
          EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
         
          coty = 1.0/tan(y); // if GridRank = 1, assume theta = pi/2

          dtxinv = dtFixed/x;
          dU[iS1][n]  += dtxinv*(2.0*p + rho*(vy*vy+vz*vz));
          dU[iS2][n]  += (GridRank > 1) ? dtxinv*(coty*p + rho*(vz*vz*coty - vx*vy)) : 0.0;
          dU[iS3][n]  += (GridRank > 2) ? -dtxinv*rho*vz*(vx+vy*coty) : 0.0;
        }
      }
    }
  }

  if ((SelfGravity) || ExternalGravity || (PointSourceGravity > 0)) {
    int igrid;
    float rho, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  //	  fprintf(stderr, "glad youre calling me");	  
	  gx = AccelerationField[0][igrid];
	  gy = (GridRank > 1) ? (AccelerationField[1][igrid]) : 0;
	  gz = (GridRank > 2) ? (AccelerationField[2][igrid]) : 0;

	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];

	  dU[iS1  ][n] += dtFixed*gx*rho;
	  dU[iS2  ][n] += dtFixed*gy*rho;
	  dU[iS3  ][n] += dtFixed*gz*rho;
	  dU[iEtot][n] += dtFixed*rho*(gx*vx + gy*vy + gz*vz);
	}
      }
    }
  }

  /*
  if (SelfGravity && GridRank == 1) {
    // calculate gravitational field
    float pi = 4.0*atan(1.0);
    float *gr = new float[GridDimension[0]-2*NumberOfGhostZones];
    FLOAT dx = CellWidth[0][0];
    for (int i = 0; i < GridDimension[0]-2*NumberOfGhostZones; i++) {
      gr[i] = 0.0;
      for (int j = GridStartIndex[0]; j < i+NumberOfGhostZones; j++) {
	gr[i] -= BaryonField[DensNum][j]*4.0*pi*pow(CellLeftEdge[0][j]+0.5*dx,2)*dx;
      }
      gr[i] -= BaryonField[DensNum][i+NumberOfGhostZones]*2.0*pi*pow(CellLeftEdge[0][i+NumberOfGhostZones]+0.5*dx,2)*dx;
      gr[i] /= pow(CellLeftEdge[0][i+NumberOfGhostZones]+0.5*dx,2);
    }
    float rho, gx;
    float vx;
    int n=0;
    for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
      rho = BaryonField[DensNum][i];	  
      vx = BaryonField[Vel1Num][i];
      dU[iS1  ][n] += dtFixed*gr[n]*rho;
      dU[iEtot][n] += dtFixed*rho*gr[n]*vx;
    }
    delete [] gr;
  }
  */

  /* Apply external driving force */

  if (UseDrivingField) {
    int Drive1Num, Drive2Num, Drive3Num;
    if (IdentifyDrivingFields(Drive1Num, Drive2Num, Drive3Num) == FAIL) {
      printf("grid::SourceTerms: cannot identify driving fields.\n");
      return FAIL;
    }
    int igrid;
    float drivex, drivey, drivez, vx, vy, vz, rho;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  drivex = BaryonField[Drive1Num][igrid];
	  drivey = BaryonField[Drive2Num][igrid];
	  drivez = BaryonField[Drive3Num][igrid];
	  vx     = BaryonField[Vel1Num      ][igrid];
	  vy     = BaryonField[Vel2Num      ][igrid];
	  vz     = BaryonField[Vel3Num      ][igrid];
	  rho    = BaryonField[DensNum     ][igrid];

	  dU[iS1  ][n] += dtFixed*rho*drivex*DrivingEfficiency;
	  dU[iS2  ][n] += dtFixed*rho*drivey*DrivingEfficiency;
	  dU[iS3  ][n] += dtFixed*rho*drivez*DrivingEfficiency;
	  dU[iEtot][n] += dtFixed*rho*(drivex*vx + drivey*vy + drivez*vz)*DrivingEfficiency;
	}
      }
    }
  }

  /* Add centrifugal force for the shearing box */


  if ((ProblemType == 35 || ProblemType == 36 ||ProblemType == 37) && ShearingBoxProblemType !=0) {
    int igrid;
    float rho, gx, gy, gz;
    FLOAT xPos[3];
    float vels[3]; 
    int n = 0;

    int iden=FindField(Density, FieldType, NumberOfBaryonFields);
    int ivx=FindField(Velocity1, FieldType, NumberOfBaryonFields);
    int ivy=FindField(Velocity2, FieldType, NumberOfBaryonFields);
    int ivz;
    if (GridRank==3)  ivz=FindField(Velocity3, FieldType, NumberOfBaryonFields);
 
    int indexNumbers[3]={iS1,iS2,iS3};

    float A[3]={0,0,0};//Omega
    A[ShearingOtherDirection]=AngularVelocity;
    
    float lengthx=DomainRightEdge[0]-DomainLeftEdge[0]; 
    float lengthy=DomainRightEdge[1]-DomainLeftEdge[1];
    float lengthz;
    if (GridRank==3) lengthz=DomainRightEdge[2]-DomainLeftEdge[2];
    else lengthz=0.0;
    

    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {

	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[iden][igrid];
	  xPos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-lengthx/2.0;
	  xPos[1] = CellLeftEdge[1][i] + 0.5*CellWidth[1][i]-lengthy/2.0;
	  if (GridRank==3) xPos[2] = CellLeftEdge[2][i] + 0.5*CellWidth[2][i]-lengthz/2.0;
	  else xPos[2]=0;
	  
	  vels[0] = BaryonField[ivx][igrid];
	  vels[1] = BaryonField[ivy][igrid];
	  if (GridRank==3) vels[2] = BaryonField[ivz][igrid];
	  else vels[2]=0;

	  //Omega cross V

	  dU[indexNumbers[0]][n] -= dtFixed*2.0*rho*(A[1]*vels[2]-A[2]*vels[1]);
	  dU[indexNumbers[1]][n] -= dtFixed*2.0*rho*(A[2]*vels[0]-A[0]*vels[2]);
	  if (GridRank==3) dU[indexNumbers[2]][n] -= dtFixed*2.0*rho*(A[0]*vels[1]-A[1]*vels[0]);
	

	  dU[indexNumbers[ShearingBoundaryDirection]][n] += dtFixed*2.0*rho*VelocityGradient*AngularVelocity*AngularVelocity*xPos[ShearingBoundaryDirection];
	  
	  
	  dU[iEtot][n] +=  dtFixed*2.0*rho*VelocityGradient*AngularVelocity*AngularVelocity*xPos[ShearingBoundaryDirection]*vels[ShearingBoundaryDirection];
	


	  
 	}
      }
    }
  }


  return SUCCESS;
}
