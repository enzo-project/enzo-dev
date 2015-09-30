/***********************************************************************
/
/  GRID CLASS (Initialize Orszag Tang Vortex)
/
/  written by: David Collins
/  date:       2004-2013
/  modified1:
/
/  PURPOSE:  Sets up the test problem.  Results can be found in most MHD 
/            method papers, such as Collins et al 2010.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
#include "phys_constants.h"
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

  FieldType[NumberOfBaryonFields++] = Density;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if( EquationOfState == 0 )
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  if( UseMHD ){
      FieldType[NumberOfBaryonFields++] = Bfield1;
      FieldType[NumberOfBaryonFields++] = Bfield2;
      FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if ( HydroMethod == MHD_RK ){
      FieldType[NumberOfBaryonFields++] = PhiField;
  }

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  float X, Y, Vx, Vy, GasEnergy=Pressure/(Gamma-1), TotalEnergy=0; 
  int index, size=1, i,j,k, field;
  float Scale[3];

  //In order to ensure a good comparison, we use the MHD-CT initialization
  //for both MHDCT and Dedner.  Temporarily, we make this code think that MHD-CT is on.

  if ( HydroMethod == MHD_RK ){
      UseMHDCT = TRUE;
      MHD_SetupDims(); //this only sets some variables that won't be used
  }

  this->AllocateGrids();  

  for(i=0;i<GridRank;i++){
    size*=GridDimension[i];
    Scale[i]=(GridRightEdge[i]-GridLeftEdge[i])/(GridDimension[i]-2*NumberOfGhostZones);
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num);

  fprintf(stderr,"Density %"FSYM" Pressure %"FSYM" V0 %"FSYM" B0 %"FSYM" \n", DensityIn,Pressure,V0,B0);
  fprintf(stderr,"Scale: %"FSYM" %"FSYM" %"FSYM"\n", Scale[0],Scale[1],Scale[2]);
  
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

        ElectricField[field][index]=B0*( cos(4.0*pi*X)/(4.0*pi) + cos(2.0*pi*Y)/(2.0*pi) );
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

        Vx=-V0*sin(2.0*pi*Y);
        Vy=V0*sin(2.0*pi*X);

        BaryonField[DensNum][index]=DensityIn;

        if( EquationOfState == 0 ){
          TotalEnergy=GasEnergy + 0.5*DensityIn*(Vx*Vx+Vy*Vy)
            +0.5*(BaryonField[B1Num][index]*BaryonField[B1Num][index]+
                  BaryonField[B2Num][index]*BaryonField[B2Num][index]+
                  BaryonField[B3Num][index]*BaryonField[B3Num][index]);
        BaryonField[TENum][index]=TotalEnergy/DensityIn;

        }
        if( DualEnergyFormalism )
          BaryonField[GENum][index]=GasEnergy;
        BaryonField[Vel1Num][index]=Vx;
        BaryonField[Vel2Num][index]=Vy;
        BaryonField[Vel3Num][index]=0.0;

      }

  if ( HydroMethod == MHD_RK ){
      //Clean up.
      UseMHDCT = FALSE;
      for ( field=0; field<3; field++){
          delete [] MagneticField[field];
          delete [] ElectricField[field];
          MagneticField[field] = NULL;
          ElectricField[field] = NULL;
      }
  }
  
  return SUCCESS;
}
