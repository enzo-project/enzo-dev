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

int QuantumGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int grid::LightBosonInitializeGrid(float CenterPosition)
{  

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = RePsi;
  FieldType[NumberOfBaryonFields++] = ImPsi;
  FieldType[NumberOfBaryonFields++] = FDMDensity;


  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }


  int size = 1, index, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];


  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (QuantumGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  int RePsiNum, ImPsiNum, FDMDensNum;
  RePsiNum = FindField(RePsi, FieldType, NumberOfBaryonFields);
  ImPsiNum = FindField(ImPsi, FieldType, NumberOfBaryonFields);
  FDMDensNum = FindField(FDMDensity, FieldType, NumberOfBaryonFields);


  float coef = (5.9157166856e27*TimeUnits/pow(LengthUnits,2));

  FLOAT x;
  float xv;

  //float a = 1./(1.+2.);
  float alpha = 1./500;
  float initialtime=0.0;
  float sumsquare = pow(alpha,2)+pow(coef*initialtime,2);

  float a = 0.1;
  float expa, expb;
  float pi = 3.1415927;

  int i;
  for (int k = 0; k < GridDimension[2]; k++) {
  for (int j = 0; j < GridDimension[1]; j++) {
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - CenterPosition;
    xv = CellLeftEdge[0][i] - CenterPosition;

    index = GRIDINDEX_NOGHOST(i,j,k); 

    // set up an initial Gaussian density profile, v=0

// 1) Gaussian Density Field
      //BaryonField[iden ][index] = 1*exp(-alpha*pow(x,2)/sumsquare)/pow(sumsquare,0.5)/5; 
      //BaryonField[RePsiNum][index] = sqrt(BaryonField[iden ][index]);
      //BaryonField[ImPsiNum][index] = 0;
      //BaryonField[ivx ][index] = 1.e-4*exp(-alpha*pow(x,2)/sumsquare); // Gaussian Density Field
      //BaryonField[ivx  ][index] = xv/sumsquare*coef*coef*initialtime;
// 2) Fresnel solution
    /*if( x>0){
      BaryonField[iden ][index] = 6.;
      BaryonField[RePsiNum][index] = sqrt(6.);
      BaryonField[ImPsiNum][index] = 0;
    }else if (x==0.){
      BaryonField[iden ][index] = 5.5;
      BaryonField[RePsiNum][index] = sqrt(5.5);
      BaryonField[ImPsiNum][index] = 0;
    }else{
      BaryonField[iden ][index] = 5.;
      BaryonField[RePsiNum][index] = sqrt(5.);
      BaryonField[ImPsiNum][index] = 0;
    }

// 3) Zeldovich Test
      //BaryonField[iden ][index] = (10.0+5.0*cos(2*3.1415927*(x-0.5)))/10.;
      //BaryonField[iden ][index] = 1.0+2.*exp(-1.*pow((x-0.7)/0.05,2))+2.*exp(-1.*pow((x-0.3)/0.05,2));
     /*BaryonField[iden  ][index] = 1;
     if ((x<0.5)&&(x>0.4)){
          BaryonField[iden  ][index] = 1+cos(10*3.1415926536*(x-0.45));
      }else if ((x>0.5)&&(x<0.6)){
          BaryonField[iden  ][index] = 1+cos(10*3.1415926536*(x-0.55));
      }
      BaryonField[ivx  ][index] = 0;
      if ((BaryonField[iden ][index]>(1+1e-6))&&(xv<0.5)){
      	  BaryonField[ivx  ][index] = 1e0;
      }else if ((BaryonField[iden ][index]>(1+1e-6))&&(xv>0.5)){
      	  BaryonField[ivx  ][index] = -1e0;
      }*/

//  4) Two colliding Gaussian packets
    expa = exp(-alpha*pow((x+a),2)/sumsquare/2.);
    expb = exp(-alpha*pow((x-a),2)/sumsquare/2.);
    float k1 = 2*pi*10;
    float k2 = 2*pi*10;
    float rho = expa*expa + expb*expb + expa*expb*2*cos(2*k*x);

    BaryonField[RePsiNum][index] = expa*cos(k1*(x+a)) + expb*cos(k2*(x-a));
    BaryonField[ImPsiNum][index] = expa*sin(k1*(x+a)) - expb*sin(k2*(x-a));
    BaryonField[FDMDensNum][index] = BaryonField[RePsiNum][index]*BaryonField[RePsiNum][index] + BaryonField[ImPsiNum][index]*BaryonField[ImPsiNum][index];
    //BaryonField[iden][index] = BaryonField[RePsiNum][index]*BaryonField[RePsiNum][index] + BaryonField[ImPsiNum][index]*BaryonField[ImPsiNum][index];

    BaryonField[ivx][index] = coef*(k1*(expa*expa-expb*expb)-2*a*alpha/sumsquare*expa*expb*sin(2*k*x))/rho;


      //BaryonField[ivx  ][index] = 0*(xv-CenterPosition)/sumsquare*coef*coef*initialtime;
      BaryonField[ivy  ][index] = 0;
      BaryonField[ivz  ][index] = 0;
      BaryonField[ietot][index] = 100; 
  }
  }
  }

   /*if (this->ComputeQuantumPressure(Time) == FAIL) {
    ENZO_FAIL("Error in ComputeQuantumPressure!\n");
  }*/

  return SUCCESS;
}

