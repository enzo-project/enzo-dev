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

int grid::ShearingBoxStratifiedInitializeGrid(float ThermalMagneticRatio, float fraction, float ShearingGeometry, int InitialMagneticFieldConfiguration)
{

  
  /* declarations */


  
 
  int phip_num;
  NumberOfBaryonFields = 0;
  FieldType[iden=NumberOfBaryonFields++] = Density;
  FieldType[ivx=NumberOfBaryonFields++] = Velocity1;
  FieldType[ivy=NumberOfBaryonFields++] = Velocity2;
  FieldType[ivz=NumberOfBaryonFields++] = Velocity3;
 
  FieldType[ietot=NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[ieint=NumberOfBaryonFields++] = InternalEnergy;
  }
  if (UseMHD){
  FieldType[iBx=NumberOfBaryonFields++] = Bfield1;
  FieldType[iBy=NumberOfBaryonFields++] = Bfield2;
  FieldType[iBz=NumberOfBaryonFields++] = Bfield3;
  FieldType[NumberOfBaryonFields++] = PhiField;
  }
  
  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }
  

  
 

  
//  int iBx, iBy, iBz;
  int iB[3]={-1,-1,-1};
  if (UseMHD){
    iB[0]=iBx;
    iB[1]=iBy;
    if (GridRank==3){ 
      iB[2]=iBz;
    }
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  this->AllocateGrids();

  srand(110182*ProcessorNumber);


  int i,j,k; 
  int n=0; 
  float eint, v2, vx,vy,vz;

  float rhou, lenu, tempu, tu, velu;
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
  float PressureUnits = rhou*pow(velu,2);

  
  float magnitude=AngularVelocity*fraction;
  float rho=1.0;

  float c_s=1e-3;
 
  float lengthx=DomainRightEdge[0]-DomainLeftEdge[0]; 
  float lengthy=DomainRightEdge[1]-DomainLeftEdge[1];
  float lengthz=DomainRightEdge[2]-DomainLeftEdge[2];

  float h, cs, dpdrho, dpde, H, pressure;  	
  float bunit=sqrt(4.0*3.14159*rhou*velu*velu);

  FLOAT x,y,z;


  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

	FLOAT xPos[3] = {CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-lengthx/2.0,
		       CellLeftEdge[1][j] + 0.5*CellWidth[1][j]-lengthy/2.0,
		       0.0};

	
	x=xPos[0]; y=xPos[1]; z=0.0;

	if (GridRank==3) {
	  xPos[2]=CellLeftEdge[2][k] + 0.5*CellWidth[2][k]-lengthz/2.0;
	  z=xPos[2];
	}
	

	float xVel[3] ={0,0,0};

	if (ShearingBoxProblemType == 0){
	  if (x*x+y*y+z*z<0.25*ShearingGeometry*ShearingGeometry){
	    BaryonField[iden ][n]=500.0*rho;
	    
	  }
	  else if (x*x+y*y+z*z<ShearingGeometry*ShearingGeometry){ 
	    BaryonField[iden ][n]=50.0*rho;
	    
	  }
	  else BaryonField[iden ][n]=rho;
	  xVel[ShearingBoundaryDirection]=10*AngularVelocity;
	 
	}
	else if (ShearingBoxProblemType == 1){ 
	  xVel[ShearingBoundaryDirection]=magnitude*sin(xPos[ShearingOtherDirection]*2.0*ShearingGeometry*3.14156);
	  xVel[ShearingVelocityDirection]=magnitude/3.*sin(xPos[ShearingVelocityDirection]*2.0*ShearingGeometry*3.14156);
	  BaryonField[iden ][n] = rho;

	  
	}

	float rhoActual=BaryonField[iden ][n];
	pressure=c_s*c_s*rhoActual/Gamma;
	float realpressure=pressure*PressureUnits;  
	float InitialBField=sqrt((8*3.14159*realpressure/(ThermalMagneticRatio)))/bunit;

	eint=0.0;
	
	if (HydroMethod == MHD_RK || HydroMethod == HD_RK) 
	  EOS(pressure, rhoActual, eint, h, cs, dpdrho, dpde, EOSType, 1);
	else eint=pressure/(rhoActual*(Gamma-1.0));

	xVel[ShearingVelocityDirection]+=(-(xPos[ShearingBoundaryDirection])*AngularVelocity*VelocityGradient);
        
	BaryonField[ivx  ][n] = xVel[0];
	BaryonField[ivy  ][n] = xVel[1];
	BaryonField[ivz  ][n] = xVel[2];

	
	
	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}

	v2=xVel[0]*xVel[0]+xVel[1]*xVel[1]+xVel[2]*xVel[2];
	
	if (UseMHD){
	  BaryonField[iBx][n] = 0.0;
	  BaryonField[iBy][n] = 0.0;
	  BaryonField[iBz][n] = 0.0;
	  
	  if (InitialMagneticFieldConfiguration == 0) 
	    BaryonField[iB[ShearingOtherDirection]][n] = InitialBField;
	  else if (InitialMagneticFieldConfiguration == 1) 
	    BaryonField[iB[ShearingOtherDirection]][n] = InitialBField*sin(2*3.14159*xPos[ShearingBoundaryDirection]);
	  
	  BaryonField[ietot][n] = eint + 0.5*v2
	    + 0.5*(BaryonField[iB[ShearingOtherDirection]][n]*
		   BaryonField[iB[ShearingOtherDirection]][n])/rhoActual;
	}
	else BaryonField[ietot][n] = eint + 0.5*v2;
	
	  
      } // end loop over grid
    }
  }
 

 
  return SUCCESS;

 

 
}
