/***********************************************************************
/
/  INITIALIZE A SHEARING BOX TEST ON A GRID
/
/  written by: Fen Zhao
/  date:       June, 2009
/  modified1:
/
/  PURPOSE:
/    Set up exither an advecting sphere or the standard shearing box simluation
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "hydro_rk/EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double Gaussian(double cs);

int grid::ShearingBox2DInitializeGrid(float ThermalMagneticRatio, float fraction, float ShearingGeometry, int InitialMagneticFieldConfiguration)
{

 // 1 to 5 :2D shearing box
    // 1 = hydro shearing box
    // 2 = KH instablity (v_y= Omega*sin(2.*3.1415927*ShearingGeometry*y)*fraction)
    // 3 = sheared sphere (radius= ShearingGeometry)
    // 4 = rotating sphere, uniform density  (radius= ShearingGeometry)
    // 5 = rotating sphere, ramped density  (radius= ShearingGeometry)


  float Pi= 3.14159;



 /* declarations */

  int phip_num;
  NumberOfBaryonFields = 0;
  FieldType[iden=NumberOfBaryonFields++] = Density;
  FieldType[ivx=NumberOfBaryonFields++] = Velocity1;
  FieldType[ivy=NumberOfBaryonFields++] = Velocity2;
 
  FieldType[ietot=NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[ieint=NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[iBx=NumberOfBaryonFields++] = Bfield1;
    FieldType[iBy=NumberOfBaryonFields++] = Bfield2;
  
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  
  if(UseDivergenceCleaning) {
    FieldType[NumberOfBaryonFields++] = Phi_pField;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }


  this->AllocateGrids();

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  float MagneticUnits = sqrt(4.0*Pi*DensityUnits)*VelocityUnits;
  
  /* Problem parameters */
  float rho = 1.0;
  float cs = 1e-3;
  float pres = rho*cs*cs;
  float Bnaught = 0.0;
  float c_s=1e-3;

  float rhou, lenu, tempu, tu, velu;
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
  float PressureUnits = rhou*pow(velu,2);
  float bunit=sqrt(4.0*3.14159*rhou*velu*velu);
 
  const float q = VelocityGradient;
  const float Omega = AngularVelocity;

 


  const FLOAT Lx = DomainRightEdge[0] - DomainLeftEdge[0];
  const FLOAT Ly = DomainRightEdge[1] - DomainLeftEdge[1];
  /* Set up the background */

  int n = 0;  
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++, n++) {

	FLOAT x;
	if (ShearingBoxProblemType==1)
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	else
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-Lx/2;

	FLOAT y;
	if (ShearingBoxProblemType==3 || ShearingBoxProblemType==4 || 
	    ShearingBoxProblemType==5)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -Ly/2;
	else
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];


	
	BaryonField[iden][n] = rho;


	float vx = 0.0;
	float vy = -x*q*Omega;
	float vz = 0.0;

	BaryonField[ivx][n] = vx;
	BaryonField[ivy][n] = vy;
	BaryonField[ivz][n] = vz;



	if (ShearingBoxProblemType==2){  //KH instability
	  float deltaRho=1.0;
	  FLOAT x1 = 0.3*Lx/2.0;
	  FLOAT delx = x1 * 0.1;
  	  FLOAT ramp =  1./((1.+exp(-2/delx*(x+x1)))*(1.+exp(-2/delx*(x1-x))));
	  float drho=deltaRho*ramp;
	 
	  BaryonField[ivy][n]+=(x*q*Omega)*ramp;
	  BaryonField[ivx][n]+=Omega*sin(2.*3.1415927*ShearingGeometry*y)*fraction;
	  BaryonField[iden][n]+= drho;
 
	}


	if (ShearingBoxProblemType==3){  //sheared sphere
	  float radius= ShearingGeometry;
	  float deltaRho=1.0;
	  FLOAT delx = radius * 0.1;
	  float r=sqrt(x*x+y*y);
	  FLOAT ramp =  1./(1.+exp(-2/delx*(radius-r)));
	  float drho=deltaRho*ramp;
	  
	  BaryonField[iden][n] += drho;
	 
	}


	if (ShearingBoxProblemType==4){  //rotating sphere, no density difference
	  float radius= ShearingGeometry;
	  //float deltaRho=1.0;
	  FLOAT delx = radius * 0.1;
	  float r=sqrt(x*x+y*y);
	  float theta=atan2(y,x);
	  //if (theta<0) theta=theta+3.1415927*2;

	  float vx=r*sin(theta);
	  float vy=r*cos(theta);

	  FLOAT ramp =  1./(1.+exp(-2/delx*(radius-r)));
	  float dvx=vx*ramp;	  
	  float dvy=vy*ramp;
	  
	  BaryonField[ivx][n] += dvx;
	  BaryonField[ivy][n] += dvy;
	 

	}

	if (ShearingBoxProblemType==5){  //rotating sphere, density difference
	  FLOAT radius=(FLOAT) ShearingGeometry;
	  float deltaRho=1.0;
	  FLOAT delx = radius * 0.1;
	
	  float r=sqrt(x*x+y*y);
	  float theta=atan(y/x);

	  float vx=r*sin(theta);
	  float vy=r*cos(theta);
	  
	  FLOAT ramp =  1./(1.+exp(-2/delx*(radius-r)));
	  float dvx=vx*ramp;	  
	  float dvy=vy*ramp;
	  float drho=deltaRho*ramp;

	  BaryonField[ivx][n] += dvx;
	  BaryonField[ivy][n] += dvy;
	  BaryonField[iden][n] += drho;
	 

	}


	vx= BaryonField[ivx][n];
	vy =BaryonField[ivy][n];
	vz= BaryonField[ivz][n];

	float eint, h, cs_temp, dpdrho, dpde;
	EOS(pres, BaryonField[iden][n], eint, h, cs_temp, dpdrho, dpde, EOSType, 1);
  	//EOS(pres, rho, eint, h, cs_temp, dpdrho, dpde, EOSType, 1);

	BaryonField[ietot][n] = eint + 0.5*(vx*vx + vy*vy + vz*vz);

	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}

	if (HydroMethod == MHD_RK) {

	  float rhoActual=BaryonField[iden ][n];
	  float pressure=c_s*c_s*rhoActual/Gamma;
	  float realpressure=pressure*PressureUnits;  
	  float InitialBField=sqrt((8*3.14159*realpressure/(ThermalMagneticRatio)))/bunit;
	  if (InitialMagneticFieldConfiguration == 0) Bnaught = InitialBField;
	  else if (InitialMagneticFieldConfiguration == 1) Bnaught = InitialBField*sin(2*3.14159*x);


	  BaryonField[iBz][n] = Bnaught;
	  BaryonField[ietot][n] += 0.5 * pow(Bnaught,2) / rho;
	}	

      } // end loop over grid
    }
  }

 

    
  if (ShearingBoxProblemType == 1){
    /* ProblemType 1: Vortex wave.
       Reference: B. M. Johnson & C. F. Gammie, ApJ, 2005, 626, 978. */

    const FLOAT kx0 = (-8.0*2.0*Pi/Lx)/(8.0);
    const FLOAT ky = 2.0*2.0*Pi/Ly;
    const float vx0 = 1e-4; // in unit of cs
    n = 0;  
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++, n++) {
	  int igrid = i + (j + k * GridDimension[1]) * GridDimension[0];
	  FLOAT x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  FLOAT y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  
	  float dvx = vx0 * cs * cos(kx0 * x + ky * y);
	  float dvy = -kx0 / ky * dvx;
	  float drho = rho/(cs*ky) * (-2.0*kx0*kx0*q*Omega + 2.0*(q-1.0)*Omega) * vx0 * sin(kx0*x + ky*y);
	  
	  //BaryonField[iden][igrid] += drho;
	 
	  BaryonField[ivx ][igrid] += dvx;
	  BaryonField[ivy ][igrid] += dvy;	  

	}
      }}}
 
 
 
    //PrintToScreenBoundaries(BaryonField[ivy], "Vy (Initial)\n" ,2, 0);
    // PrintToScreenBoundaries(BaryonField[iden], "Dens (Initial)\n",  2, 0);
    

  return SUCCESS;
 
}
