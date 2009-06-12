#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double Gaussian(double cs);

int grid::MRICollapseInitializeGrid(float AngularVelocity, float VelocityGradient, float ThermalMagneticRatio, float fraction, float radius)
{
  /* declarations */


  int phip_num;
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[ietot=NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[ieint=NumberOfBaryonFields++] = InternalEnergy;
  }
  FieldType[NumberOfBaryonFields++] = Bfield1;
  FieldType[NumberOfBaryonFields++] = Bfield2;
  FieldType[NumberOfBaryonFields++] = Bfield3;
  FieldType[NumberOfBaryonFields++] = PhiField;
  
  if(UseDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
    FieldType[NumberOfBaryonFields++] = DebugField;  
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
    }
  }


  srand(110182);
  float magnitude=AngularVelocity*fraction;

  int i,j,k; int n=0; float eint, v2, vx,vy,vz;

  float rhou, lenu, tempu, tu, velu;
   GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
  float PressureUnits = rhou*pow(velu,2);

  

  float rho=1.0;
  
  float realRho=rho*rhou;
  float realRadius=radius*lenu;
  float Mass=realRho*(4*3.14159/3)*pow(realRadius,3);
 
  float G=6.67259e-8;

  float realPressurec=(3.0/(8.0*3.14159))*G*Mass*Mass/pow(realRadius,4);
  float pressurec=realPressurec/PressureUnits;



  fprintf(stdout, "initial values %g %g %g %g  (%g %g %g)%g %g\n", realRho, realRadius, Mass, G, (3.0/(8.0*3.14159))*G,Mass*Mass,pow(realRadius,4), realPressurec, pressurec);

  
  float h, cs, dpdrho, dpde, H, pressure;


  float bunit=sqrt(4.0*3.14159*rhou*velu*velu);
  float InitialBField=sqrt((8*3.14159*realPressurec/(ThermalMagneticRatio*bunit)));
  //InitialBField=1.121e-7

  //float InitialBField=sqrt((2*pressure/(ThermalMagneticRatio)));
  


  float dsize[3]={DomainRightEdge[0]-DomainLeftEdge[0],
		 DomainRightEdge[1]-DomainLeftEdge[1],
		 DomainRightEdge[2]-DomainLeftEdge[2]};
 
  float omega;
 
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

	FLOAT x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-0.5*dsize[0];
	FLOAT y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j]-0.5*dsize[1];
	FLOAT z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k]-0.5*dsize[2];
	
	FLOAT r=pow(x*x+y*y+z*z, 0.5);

	pressure=pressurec*(1-pow(r/radius,2));

	if (r<=radius){
	  BaryonField[iden ][n] = rho;
	  pressure=pressurec*(1-pow(r/radius,2));
	}
	else{
	  //BaryonField[iden ][n] = 0;
	  BaryonField[iden ][n] = rho*tiny_number;
	  pressure=pressurec*tiny_number*tiny_number;
	}
	eint=0.0;
	if (DualEnergyFormalism) {
	  EOS(pressure, BaryonField[iden ][n], eint, h, cs, dpdrho, dpde, EOSType, 1);
	}
  
	if (r<=radius){ 
	  vx=magnitude*(rand()-RAND_MAX/2)/RAND_MAX;
	  vy=magnitude*(rand()-RAND_MAX/2)/RAND_MAX; 
	  vz=magnitude*(rand()-RAND_MAX/2)/RAND_MAX;
	  

	  vx=vx - AngularVelocity*(y);
	  vy=vy + AngularVelocity*(x);
	  vz=0;
	  
	}
 	else{
 	  vx=0.0;
 	  vy=0.0;
 	  vz=0.0;
 	}


	BaryonField[ivx  ][n] = vx;
	BaryonField[ivy  ][n] = vy;
	BaryonField[ivz  ][n] = vz;
	

	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}

	//H=cs/AngularVelocity;
	
	BaryonField[iBx][n] = 0.0;
	//BaryonField[iBx][n] = InitialBField;
	BaryonField[iBy][n] = 0.0;
	//BaryonField[iBz][n] = InitialBField;
	BaryonField[iBz][n] = 0.0;
	BaryonField[iPhi][n] = 0.0;
	if(UseDivergenceCleaning) BaryonField[phip_num][n] = 0.0;




	v2=vx*vx + vy*vy + vz*vz;

	BaryonField[ietot][n] = eint + 0.5*v2
	  + 0.5*(BaryonField[iBz][n]*BaryonField[iBz][n])/rho;
	



      } // end loop over grid
    }
  }



//   SetShearingBoxExternalBoundaries();
//  n = getIndex(4,4,4);
//    fprintf(stdout, "values %f %f %f\n", BaryonField[iden ][n], BaryonField[ieint][n],  	BaryonField[ietot][n]);

  return SUCCESS;
}





