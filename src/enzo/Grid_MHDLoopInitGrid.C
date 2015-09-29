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
#include "DebugTools.h"
#include "phys_constants.h"


int grid::MHDLoopInitGrid(float LoopDensity,float Pressure, float Vx, float Vy, float Vz, float B0, FLOAT R0, 
                          FLOAT Center[], int CurrentAxis){ 


  fprintf(stderr,"GridDim %"ISYM" %"ISYM" %"ISYM"\n",GridDimension[0],GridDimension[1],GridDimension[2]);
  int field=0;

  FieldType[NumberOfBaryonFields++] = Density;
  if( EquationOfState == 0 ) FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  if( UseMHD ){
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if( HydroMethod == MHD_RK ){
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  if(DualEnergyFormalism) FieldType[NumberOfBaryonFields++] = InternalEnergy;

  int Eeng, Eden, Ev[3], Egas, BxNum = 0, ByNum = 1, BzNum = 2;
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], Ev[2], Eeng,
        BxNum, ByNum, BzNum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }


  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"======== Loop ===================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");
  fprintf(stderr,"=================================\n");



  float Pi = 3.14159265, One=1.0;
  float R,X, Y, GasEnergy=Pressure/(Gamma-1), LoopTotalEnergy=0; 
  int index, size=1, i,j,k, Three=3,TENum=1; 
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
  for( i=GridRank; i<3; i++){
    Scale[i] = 0;
    Center[i] = 0;
  }

  //Vector Potential. 
  //Due to the similarity in centering, and lack of foresight in naming,
  //I'm using the Electric Field as a Vector Potential to initialize the Magnetic Field.  

  for(field=0;field<3;field++)
    for( k=0;k<ElectricDims[field][2];k++)
      for( j=0;j<ElectricDims[field][1];j++)
        for( i=0;i<ElectricDims[field][0];i++){

          index=i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);

          if( field== CurrentAxis ){

            switch( CurrentAxis ){
              case 0:
                X=(j-GridStartIndex[1])*Scale[1]-(Center[1]-DomainLeftEdge[1]);
                Y=(k-GridStartIndex[2])*Scale[2]-(Center[2]-DomainLeftEdge[2]);
                break;
              case 1:
                X=(k-GridStartIndex[2])*Scale[2]-(Center[2]-DomainLeftEdge[2]);
                Y=(i-GridStartIndex[0])*Scale[0]-(Center[0]-DomainLeftEdge[0]);
                break;
              case 2:
                X=(i-GridStartIndex[0])*Scale[0]-(Center[0]-DomainLeftEdge[0]);
                Y=(j-GridStartIndex[1])*Scale[1]-(Center[1]-DomainLeftEdge[1]);
                break;
            }//Current Axis switch.
            R=sqrt(X*X+Y*Y);
            ElectricField[field][index]=(R<R0)? B0*(R0-R):0.0;

          }else{
            ElectricField[field][index]=0.0;

          }

    }
  

  //
  //Curl is B=Curl(A)
  //


  //the last argument indicates that this isn't a time update.  See the source for details.
  if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL )
    {fprintf(stderr," error occored in MHD_Curl\n"); return FAIL;}


  if( this->CenterMagneticField() == FAIL ) 
    {fprintf(stderr," error with CenterMagneticField , second call\n");return FAIL;}

  //Set the rest of the fields.
  for(k=0;k<GridDimension[2];k++)
    for(j=0;j<GridDimension[1];j++)
      for(i=0;i<GridDimension[0];i++){

        index=i+GridDimension[0]*(j+GridDimension[1]*k);
        X=(i-GridStartIndex[0])*Scale[0];
        Y=(j-GridStartIndex[1])*Scale[1];
        
        LoopTotalEnergy=GasEnergy + 0.5*(Vx*Vx + Vy*Vy + Vz*Vz)
          +0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index]+
                BaryonField[ByNum][index]*BaryonField[ByNum][index]+
                BaryonField[BzNum][index]*BaryonField[BzNum][index])/LoopDensity;

        BaryonField[Eden][index]=LoopDensity;
        if( EquationOfState == 0 ) BaryonField[Eeng][index]=LoopTotalEnergy;
        if (DualEnergyFormalism)
          BaryonField[Egas][index]=GasEnergy/LoopDensity;
        BaryonField[ Ev[0] ][index]=Vx;
        BaryonField[ Ev[1] ][index]=Vy;
        BaryonField[ Ev[2] ][index]=Vz;

      }
  
  if(HydroMethod == MHD_RK){  
      UseMHDCT = FALSE;
    for(field=0;field<3;field++){

      delete MagneticField[field];
      MagneticField[field] = NULL;
      delete ElectricField[field];
      ElectricField[field] = NULL;
    }
    for(int field=0; field<3; field++){
      MagneticSize[field] = 0;
      ElectricSize[field] = 0;
    }
  }//HydroRK

  return SUCCESS;

}
