#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
//start

extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
					float *ex, float *ey, float *ez, 
					float *dx, float *dy, float *dz, 
					int *idim, int *jdim, int *kdim,
					int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
					float *dt, int *method);


int grid::MHDOrszagTangInitGrid(float DensityIn,float Pressure, float V0, float B0 ){ 


  //Every processor needs to know this for every grid,
  //WHETHER OR NOT IT HAS THE DATA.
  

  NumberOfBaryonFields = 0;
  fprintf(stderr,"GridDim %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
  FieldType[NumberOfBaryonFields++] = Density;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if( EquationOfState == 0 )
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"======== tang ===================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");

  float Pi = 3.14159265, One=1.0;
  float X, Y, Vx, Vy, GasEnergy=Pressure/(Gamma-1), TotalEnergy=0; 
  int index, size=1, i,j,k, field;
  float Scale[3];
  this->AllocateGrids();  

  for(i=0;i<GridRank;i++){
    size*=GridDimension[i];
    Scale[i]=(GridRightEdge[i]-GridLeftEdge[i])/(GridDimension[i]-2*NumberOfGhostZones);

  }
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  fprintf(stderr,"Density %f Pressure %f V0 %f B0 %f \n", DensityIn,Pressure,V0,B0);
  fprintf(stderr,"Scale: %f %f %f\n", Scale[0],Scale[1],Scale[2]);
  

  //Vector Potential. 
  //Due to the similarity in centering, and lack of foresight in naming,
  //I'm using the Electric Field as a Vector Potential to initialize the Magnetic Field.  

  field=2;
  for( k=0;k<ElectricDims[field][2];k++)
    for( j=0;j<ElectricDims[field][1];j++)
      for( i=0;i<ElectricDims[field][0];i++){
	index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
	X=(i-GridStartIndex[0])*Scale[0];
	Y=(j-GridStartIndex[1])*Scale[1];

	ElectricField[field][index]=B0*( cos(4*Pi*X)/(4*Pi) + cos(2*Pi*Y)/(2*Pi) );
      }
  

  //
  //Curl is B=Curl(A)
  //


  //the last argument indicates that this isn't a time update.
  if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL )
    {fprintf(stderr," error occored in MHD_Curl\n"); return FAIL;}

  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField\n");
    return FAIL;
  }



  for(k=0;k<GridDimension[2];k++)
    for(j=0;j<GridDimension[1];j++)
      for(i=0;i<GridDimension[0];i++){
	index = i+GridDimension[0]*(j+GridDimension[1]*k);
	X=(i-GridStartIndex[0]+0.5)*Scale[0];
	Y=(j-GridStartIndex[1]+0.5)*Scale[1];

	Vx=-V0*sin(2*Pi*Y);
	Vy=V0*sin(2*Pi*X);

	BaryonField[DensNum][index]=DensityIn;

	if( EquationOfState == 0 ){
	  TotalEnergy=GasEnergy + 0.5*DensityIn*(Vx*Vx+Vy*Vy)
	    +0.5*(CenteredB[0][index]*CenteredB[0][index]+
		  CenteredB[1][index]*CenteredB[1][index]+
		  CenteredB[2][index]*CenteredB[2][index]);
	BaryonField[TENum][index]=TotalEnergy/DensityIn;

	}
	if( DualEnergyFormalism )
	  BaryonField[GENum][index]=GasEnergy;
	BaryonField[Vel1Num][index]=Vx;
	BaryonField[Vel2Num][index]=Vy;
	BaryonField[Vel3Num][index]=0.0;

      }
  

  return SUCCESS;
}
