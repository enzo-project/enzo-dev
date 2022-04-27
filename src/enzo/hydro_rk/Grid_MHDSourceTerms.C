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
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"
#include "hydro_rk/SuperNova.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int field, int farray[], int numfields);

void mt_init(unsigned_int seed);
unsigned_long_int mt_random();

int grid::MHDSourceTerms(float **dU, float min_coeff)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, PhiNum, CRNum;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
				   TENum, B1Num, B2Num, B3Num, PhiNum, CRNum);
  if (CRModel) {
    if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
      ENZO_FAIL("Cannot Find Cosmic Rays");
  }

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

  if (CRModel){

    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    float coef = 0.5;
    float dtdx = coef*dtFixed/CellWidth[0][0]/a;
    float dtdy = (GridRank > 1) ? coef*dtFixed/CellWidth[1][0]/a : 0.0;
    float dtdz = (GridRank > 2) ? coef*dtFixed/CellWidth[2][0]/a : 0.0;

    FLOAT divVdt, dtdEcrdx, dtdEcrdy, dtdEcrdz, dHeatCR, rho, inv_sqrt_rho;
    FLOAT va_x, va_y, va_z;
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

          // CR energy density gradient                                                                                            
          dtdEcrdx = dtdx * (BaryonField[CRNum][ip1]-BaryonField[CRNum][im1]);
          dtdEcrdy = dtdy * (BaryonField[CRNum][jp1]-BaryonField[CRNum][jm1]);
          dtdEcrdz = dtdz * (BaryonField[CRNum][kp1]-BaryonField[CRNum][km1]);

          divVdt = dtdx*(BaryonField[Vel1Num][ip1] - BaryonField[Vel1Num][im1])  +
                   dtdy*(BaryonField[Vel2Num][jp1] - BaryonField[Vel2Num][jm1]) +
                   dtdz*(BaryonField[Vel3Num][kp1] - BaryonField[Vel3Num][km1]);

          rho = BaryonField[DensNum][igrid];

	  // CR Advection and PdV work on thermal gas
	  dU[iCR][n] += (CRgamma - 1.0) * (BaryonField[Vel1Num][igrid] * dtdEcrdx + 
		     BaryonField[Vel2Num][igrid]*dtdEcrdy + BaryonField[Vel3Num][igrid] * dtdEcrdz);
	  dU[iEtot][n] += (CRgamma - 1.0) * BaryonField[CRNum][igrid] * divVdt; 
	  if (DualEnergyFormalism)
	    dU[iEint][n] += (CRgamma - 1.0) * BaryonField[CRNum][igrid] *divVdt;


	  if (CRHeating){
	    // components of Alfven velocity
            inv_sqrt_rho = 1.0 / sqrt(rho);
	    va_x = BaryonField[B1Num][igrid] * inv_sqrt_rho;
	    va_y = BaryonField[B2Num][igrid] * inv_sqrt_rho;
	    va_z = BaryonField[B3Num][igrid] * inv_sqrt_rho;

	    // Calculating the heating rate of cosmic rays on the thermal gas:                                                       
	    dHeatCR = (CRgamma - 1.0)*(va_x*dtdEcrdx + va_y*dtdEcrdy +va_z*dtdEcrdz);

	    // Streaming faster than the Alfven velocity doesn't heat.
	    if (CRStreamVelocityFactor < 1)
	      dHeatCR *= CRStreamVelocityFactor;

	    dU[iCR][n]   -= fabs(dHeatCR);
	    if (DualEnergyFormalism)
	        dU[iEint][n] += fabs(dHeatCR);
	    dU[iEtot][n] += fabs(dHeatCR);
	  }
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

  if (UseSGSModel) {
    // if an explicit filtering operation should be used, otherwise
    // grid-scale quantities are used
    if (SGSFilterWidth > 1.) {
      if (this->SGSUtil_FilterFields() == FAIL) {
        fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_FilterFields.\n");
        return FAIL;
      }

      // if the partial derivatives of primitive variables are required
      // in the calculation of the SGS models
      if (SGSNeedJacobians) {
        // velocity Jacobian
        if (this->SGSUtil_ComputeJacobian(JacVel,FilteredFields[1],FilteredFields[2],FilteredFields[3]) == FAIL) {
          fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_ComputeJacobian(Vel).\n");
          return FAIL;
        }
        // magnetic field Jacobian
        if (this->SGSUtil_ComputeJacobian(JacB,FilteredFields[4],FilteredFields[5],FilteredFields[6]) == FAIL) {
          fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_ComputeJacobian(B).\n");
          return FAIL;
        }
      }

      // Scale-similarity type models need filtered mixed terms, such as flt(u_i B_j), etc.
      if (SGSNeedMixedFilteredQuantities) {
        if (this->SGSUtil_ComputeMixedFilteredQuantities() == FAIL) {
          fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_ComputeMixedFilteredQuantities().\n");
          return FAIL;
        }
      }

    } else {
      /* we don't need a special check for SGSNeedJacobians here as all models apart
       * from the scale-similarity model need Jacbobians and the scale-similarity model
       * always has SGSFilterWidth > 1.
       */
      if (this->SGSUtil_ComputeJacobian(JacVel,BaryonField[Vel1Num],BaryonField[Vel2Num],BaryonField[Vel3Num]) == FAIL) {
        fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_ComputeJacobian(Vel).\n");
        return FAIL;
      }
      if (this->SGSUtil_ComputeJacobian(JacB,BaryonField[B1Num],BaryonField[B2Num],BaryonField[B3Num]) == FAIL) {
        fprintf(stderr, "grid::MHDSourceTerms: Error in SGSUtil_ComputeJacobian(B).\n");
        return FAIL;
      }
    }

    if (this->SGS_AddMomentumTerms(dU) == FAIL) {
      fprintf(stderr, "grid::MHDSourceTerms: Error in SGS_AddMomentumTerms(dU).\n");
      return FAIL;
    }

    if (this->SGS_AddEMFTerms(dU) == FAIL) {
      fprintf(stderr, "grid::MHDSourceTerms: Error in SGS_AddEMFTerms(dU).\n");
      return FAIL;
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
  
  if(UseMagneticSupernovaFeedback) {
 
    int n, active_x, active_y, center_i, center_j, center_k, num_sn_cells_x, num_sn_cells_y, num_sn_cells_z; 
    snsf_source_terms S;
    float dx, dy, dz, dist_to_sn, magnetic_energy_density;
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;

    if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, Time) == FAIL){
      fprintf(stderr, "Error in GetUnits.\n");
      return FAIL;
    }

    // Converting radius from parsecs to cm, then internal units   
    float sn_radius = MagneticSupernovaRadius * pc_cm / LengthUnits;

    active_x = GridDimension[0] - 2*NumberOfGhostZones; 
    active_y = GridDimension[1] - 2*NumberOfGhostZones;

    // the number of additional cells that receive magnetic supernova feedback. 
    // assuming the supernova is in the center of one cell
    num_sn_cells_x = (int) ((sn_radius - 0.5*CellWidth[0][0])/ CellWidth[0][0]); 
    num_sn_cells_y = (int) ((sn_radius - 0.5*CellWidth[1][0])/ CellWidth[1][0]);
    num_sn_cells_z = (int) ((sn_radius - 0.5*CellWidth[2][0])/ CellWidth[2][0]);
			    
    for (std::vector<SuperNova>::iterator current_sn = this->MagneticSupernovaList.begin(); 
           current_sn != this->MagneticSupernovaList.end(); current_sn++){

      // find index of the cell nearest to the supernova center
      // assuming that supernova in the center of that cell
      center_i  = (int)((current_sn->getPosition()[0] - GridLeftEdge[0]) / CellWidth[0][0]);  
      center_j  = (int)((current_sn->getPosition()[1] - GridLeftEdge[1]) / CellWidth[1][0]);
      center_k  = (int)((current_sn->getPosition()[2] - GridLeftEdge[2]) / CellWidth[2][0]);

      for(int k = center_k - num_sn_cells_z; k <= center_k + num_sn_cells_z; k++){
	for(int j = center_j - num_sn_cells_y; j <= center_j + num_sn_cells_y; j++){
	  for(int i = center_i - num_sn_cells_x; i <= center_i + num_sn_cells_x; i++){
	    
	    // only add magnetic feedback on the active grid cells
	    if ((k >= GridStartIndex[2]) && (k <= GridEndIndex[2]) && 
		(j >= GridStartIndex[1]) && (j <= GridEndIndex[1]) &&
		(i >= GridStartIndex[0]) && (i <= GridEndIndex[0])){
	      
	      dx = CellWidth[0][0] * (float)(i-center_i); 
	      dy = CellWidth[1][0] * (float)(j-center_j);
	      dz = CellWidth[2][0] * (float)(k-center_k);
	 
	      dist_to_sn = sqrt(dx*dx + dy*dy + dz*dz);
	      S = current_sn->getSourceTerms(dx, dy, dz, Time);
	    
	      // solving for index n
	      // analogous to how igrid is calculated, but taking into acount Ghost Zones
	      n = (i - GridStartIndex[0])+((j-GridStartIndex[1]) + (k-GridStartIndex[2])*active_y) * active_x;
	      
	      dU[iBx][n] += S.dbx*dtFixed;
	      dU[iBy][n] += S.dby*dtFixed;
	      dU[iBz][n] += S.dbz*dtFixed;

	      dU[iEtot][n] += S.dUtot * dtFixed;
	      
	    }

	  } // End of k for-loop     
	} // End of j for-loop    
      } // End of i for-loop  
    } // End of MagneticSupernovaList loop
  } // End of UseMagneticSupernovaFeedback scope                                                                                   

  

  return SUCCESS;
}
