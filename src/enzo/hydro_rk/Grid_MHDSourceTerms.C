/***********************************************************************
/
/  GRID CLASS (COMPUTE MHD SOURCE TERMS)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#define USE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "list.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int field, int farray[], int numfields);


int grid::MHDSourceTerms(float **dU)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
				   TENum, B1Num, B2Num, B3Num, PhiNum);


  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }
  


  if (DualEnergyFormalism) {
    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    FLOAT coef = 0.5;
    FLOAT dtdx = coef*dtFixed/CellWidth[0][0]/a,
      dtdy = (GridRank > 1) ? coef*dtFixed/CellWidth[1][0]/a : 0.0,
      dtdz = (GridRank > 2) ? coef*dtFixed/CellWidth[2][0]/a : 0.0;
    float min_coeff = 0.0;
    if (UseMinimumPressureSupport) {
      min_coeff = MinimumPressureSupportParameter*
	0.32*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));
    }
    float rho, eint, p, divVdt, h, cs, dpdrho, dpde;
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

  if (Coordinate == Cylindrical) {
    float rho, etot, eint, vx, vy, vz, v2, e, h, cs, p, 
      dpdrho, dpde, coty, Bx, By, Bz, B2;
    FLOAT x, dtxinv;
    int n = 0, igrid;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = i+(j+k*GridDimension[1])*GridDimension[0];
          x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

          rho = BaryonField[DensNum][igrid];
          vx  = BaryonField[Vel1Num][igrid];
          vy  = BaryonField[Vel2Num][igrid];
          vz  = BaryonField[Vel3Num][igrid];
	  Bx  = BaryonField[B1Num  ][igrid];
	  By  = BaryonField[B2Num  ][igrid];
	  Bz  = BaryonField[B3Num  ][igrid];
	  if (DualEnergyFormalism) {
	    eint = BaryonField[ieint][igrid];
	  }
	  else {
	    etot = BaryonField[TENum][igrid];
	    v2 = vx*vx + vy*vy + vz*vz;
	    B2 = Bx*Bx + By*By + Bz*Bz;
	    eint = etot - 0.5*v2 - 0.5*B2/rho;
	  }
                  
          EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
         
          dtxinv = dtFixed/x;
          dU[iS1][n]  += dtxinv*(p + rho*vz*vz);
          dU[iS3][n]  += -dtxinv*rho*vx*vz;

	
        }
      }
    }
  }

  if (UseConstantAcceleration) {
    int igrid;
    float rho, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
	  gx = ConstantAcceleration[0];
	  gy = ConstantAcceleration[1];
	  gz = ConstantAcceleration[2];
	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];
	  
	  dU[iS1][n] += dtFixed*gx*rho;
	  dU[iS2][n] += dtFixed*gy*rho;
	  dU[iS3][n] += dtFixed*gz*rho;
	  dU[iEtot][n] += dtFixed*rho*(gx*vx + gy*vy + gz*vz);

	if (i==3 && j==3 && k==4 && GridLeftEdge[0]==0.0 && GridLeftEdge[1]==1.0)
	  printf("StermStart4 old %"GSYM" \n", dU[iS2][n])  ;
	}
      }
    }
  }

  if ((UseGasDrag != 0) && (GasDragCoefficient != 0.)) {
    int igrid;
    float rho;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
	  vx = BaryonField[Vel1Num][igrid];
	  vy = BaryonField[Vel2Num][igrid];
	  vz = BaryonField[Vel3Num][igrid];
	  
	  dU[iS1][n] -= dtFixed*GasDragCoefficient*vx*rho;
	  dU[iS2][n] -= dtFixed*GasDragCoefficient*vy*rho;
	  dU[iS3][n] -= dtFixed*GasDragCoefficient*vz*rho;
	  dU[iEtot][n] -= dtFixed*rho*(GasDragCoefficient*vx*vx + 
				       GasDragCoefficient*vy*vy + 
				       GasDragCoefficient*vz*vz);

	if (i==3 && j==3 && k==4 && GridLeftEdge[0]==0.0 && GridLeftEdge[1]==1.0)
	  printf("StermStart4 old %"GSYM" \n", dU[iS2][n])  ;
	}
      }
    }
  }


  if ((SelfGravity) || ExternalGravity || UniformGravity || PointSourceGravity) {
    int igrid;
    float rho, gx, gy, gz;
    float vx, vy, vz, vx_old, vy_old, vz_old;
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  rho = BaryonField[DensNum][igrid];
	  
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


  if (ComovingCoordinates == 1) { // add some B related cosmological expansion terms here

    int igrid;
    float rho, coef=0.;
    int n = 0;
    coef = -0.5*dadt/a;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	  //	  rho = 0.5*(BaryonField[DensNum][igrid]+OldBaryonField[DensNum][igrid]);
	  rho = BaryonField[DensNum][igrid];
	  
	  dU[iBx  ][n] += dtFixed*coef*BaryonField[B1Num][igrid];
	  dU[iBy  ][n] += dtFixed*coef*BaryonField[B2Num][igrid];
	  dU[iBz  ][n] += dtFixed*coef*BaryonField[B3Num][igrid];

	  dU[iEtot][n] -= dtFixed*coef*(BaryonField[B1Num][igrid]*BaryonField[B1Num][igrid]+
					BaryonField[B2Num][igrid]*BaryonField[B2Num][igrid]+
					BaryonField[B3Num][igrid]*BaryonField[B3Num][igrid]);

	  dU[iPhi][n] += 0.0; // Zero is the correct update, Phi is fully comoving.


	}
      }
    }
  }


  /* Apply external driving force */

  if (UseDrivingField) {

    int Drive1Num, Drive2Num, Drive3Num;
    if (IdentifyDrivingFields(Drive1Num, Drive2Num, Drive3Num) == FAIL) {
      printf("grid::SourceTerms: canot identify driving fields.\n");
      return FAIL;
    }
    int igrid;
    float drivex, drivey, drivez, vx, vy, vz, rho;
    int n = 0;
    
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i+(j+k*GridDimension[1])*GridDimension[0];

	  rho    = BaryonField[DensNum     ][igrid];
	  drivex = BaryonField[Drive1Num][igrid];
	  drivey = BaryonField[Drive2Num][igrid];
	  drivez = BaryonField[Drive3Num][igrid];
	  vx     = BaryonField[Vel1Num      ][igrid];
	  vy     = BaryonField[Vel2Num      ][igrid];
	  vz     = BaryonField[Vel3Num      ][igrid];

	  dU[iS1  ][n] += dtFixed*rho*drivex*DrivingEfficiency;
	  dU[iS2  ][n] += dtFixed*rho*drivey*DrivingEfficiency;
	  dU[iS3  ][n] += dtFixed*rho*drivez*DrivingEfficiency;
	  dU[iEtot][n] += dtFixed*rho*(drivex*vx + drivey*vy + drivez*vz +
        0.5 * dtFixed * (drivex * drivex + drivey * drivey + drivez * drivez)) * DrivingEfficiency;


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
  
if((UseSupernovaSeedFieldSourceTerms == 1)) {

  int n = 0, igrid;
  int iden=FindField(Density, FieldType, NumberOfBaryonFields);
  snsf_source_terms S;
  List<SuperNova>::Iterator *P = this->SuperNovaList.begin();
  FLOAT cell_center[3];
  FLOAT dx, dy, dz, dist_to_sn;
  int temp =1;
  int entered = 0;

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	igrid = i+(j+k*GridDimension[1])*GridDimension[0];

	cell_center[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	cell_center[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	cell_center[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	while(P != this->SuperNovaList.end()){
	  dx = P->get()->getPosition()[0] - cell_center[0];
	  dy = P->get()->getPosition()[1] - cell_center[1];
	  dz = P->get()->getPosition()[2] - cell_center[2];

	  dist_to_sn = sqrt(dx*dx + dy*dy + dz*dz);
	  if (dist_to_sn < 1.1*SupernovaSeedFieldRadius){
	    S = P->get()->getSourceTerms(dx, dy, dz, Time);
   	    double rho = BaryonField[DensNum][igrid];

	    dU[iBx][n] += S.dbx*dtFixed;
	    dU[iBy][n] += S.dby*dtFixed;
	    dU[iBz][n] += S.dbz*dtFixed;
	    dU[iEtot][n] += S.dUtot*dtFixed;

	  }
	  P = P->next();
	}// End of SuperNovaList iteration                                                                                           
      } // End of k for-loop                                                                                                         
    } // End of j for-loop                                                                                                           
  } // End of i for-loop                                                                                                             

} // End of UseSuperNovaSeedFieldSourceTerms scope                                                                                   

  

  return SUCCESS;
}
