/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR MHD Li (MHDCT) method
/
/
/  written by: dcollins
/  date:       March 19, 2013, 3:54 pm.  
/  modified1:
/
/  PURPOSE:  
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

extern "C" void FORTRAN_NAME(pde1dsolver_mhd)(float * wx, int * idim,int * nu,
           int * startindex, int * endindex, 
           float * dx, float * dtstrang,
           float * fluxBx, //remove
           float * fluxx, 
           float * fluxEx, //remove
           float * diffcoefx, 
           float * gamma, float * csmin, float * rhomin, 
           int * MHDCTDualEnergyMethod, int * MHDCTSlopeLimiter,int * RiemannSolver, 
           int * ReconstructionMethod, int * idiffusion, int * MHDCTPowellSource,
           float * tdum0, float * boxl0, float * hubb, float * zr, 
           int * nhy, int * gravityon, float * gravityx, 
           float * a, int * EquationOfState, float * SoundSpeed, int * hack2);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::SolveMHD_Li(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[])
{
    

  //Big questions: 
  //1.) Strang?  Do I want to do the full pairwise strang, or the simpler style?
  //2.) What is needed for rank < 3?
  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);
  FLOAT a[4];
  
  if(ComovingCoordinates==1)
    {
      //CosmologyComputeExpansionFactor(Time, &a[0], &a[1]) 
      CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed, &a[2], &a[3]);
    }else{
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 1.0;
    a[3] = 0.0;
  }

  int line_size = max( GridDimension[0], max(GridDimension[1], GridDimension[2]));
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int ixyz = CycleNumber % GridRank;
  int n, ii, jj, kk, index_bf, index_line;
  int hack = 0; //a flag for testing the solver.

  int nu = 6; //Remove this
  float csmin = 1e-13, rhomin = 1e-6;
  float tdum0 = 1e-13, boxl0 = 0.0, hubb=0.0, zr = 0.0; //Dummy variables.  Clean up.


  int line_width = 9;  //the number of conserved quantities.
  float * field_line = new float[line_size * line_width];
  float * flux_line = new float[line_size * line_width];
  float * colour_line = new float[line_size * NumberOfColours];
  float * gravity_line = new float[line_size];
  float * diffusion_line = new float[line_size];

  //Both of these need to be convolved with the above.
  float * flux_magnetic_line = new float[ line_size ]; //for Powell fluxes
  float * flux_def_line = new float[ line_size ]; //For dual energy.  Get rid of this.

  //DO setup strang
  //DO Comoving magnetic conversion
  float * pressure, *entropy;
  if( DualEnergyFormalism ){
    //Pressure is computed in-place for "historical reasons"
    //DO clean this up.
    pressure = BaryonField[GENum];
    for ( ii=0; ii<size; ii++)
      pressure[ii] *= BaryonField[DensNum][ii] * (Gamma - 1);
  }else{
    pressure = new float[size];
    if( EquationOfState == 0 ){
      for( ii=0; ii<size; ii++){
        pressure[ii] = (Gamma-1)*(BaryonField[TENum][ii] -
            0.5*(BaryonField[Vel1Num][ii]*BaryonField[Vel1Num][ii]+
                 BaryonField[Vel2Num][ii]*BaryonField[Vel2Num][ii]+
                 BaryonField[Vel3Num][ii]*BaryonField[Vel3Num][ii])*BaryonField[DensNum][ii]
            +0.5*(CenteredB[0][ii]*CenteredB[0][ii]+
                  CenteredB[1][ii]*CenteredB[1][ii]+
                  CenteredB[2][ii]*CenteredB[2][ii]));
      }
    }
  }
  entropy = new float[size];
  for( ii=0; ii<size; ii++){
    entropy[ii] = pressure[ii]/POW(BaryonField[DensNum][ii], Gamma - 1);
  }

  //strang loop.
  //
  int startindex, endindex;
//  for (n = ixyz; n < ixyz+GridRank; n++) {
  startindex = 3;
  endindex = GridDimension[0] - 2;
  for( n=0; n<1; n++){

    if ( 1 ){ //x loope
      for( kk = 0; kk< GridDimension[2]; kk++){
        for( jj = 0; jj<GridDimension[1]; jj++){
          for(ii = 0; ii<GridDimension[0]; ii++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ ii + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*1] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*2] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*3] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*4] = CenteredB[0][index_bf];
            field_line[ ii + line_size*5] = CenteredB[1][index_bf];
            field_line[ ii + line_size*6] = CenteredB[2][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ ii + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ ii + line_size*8] = entropy[index_bf];
            }
   
            //DO fill colour lines
            //DO fill gravity line
            gravity_line[ ii ] = 0.0;
            //DO zero flux lines
          }
            //DO diffusion term
           for( ii=1; ii<GridDimension[0]; ii++){
             diffusion_line[ii] = 0.0;
           }
            FORTRAN_NAME(pde1dsolver_mhd)(field_line, &line_size, &nu, 
                &startindex, &endindex,
                CellWidthTemp[0],  &dtFixed, 
                flux_magnetic_line, //kill this
                flux_line, 
                flux_def_line, //kill this too
                diffusion_line,
                &Gamma, &csmin, &rhomin,
                &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
                &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
                &tdum0, &boxl0, &hubb, &zr, 
                &CycleNumber, &GravityOn, gravity_line, 
                a, &EquationOfState, &IsothermalSoundSpeed, &hack);

          for(ii = 0; ii<GridDimension[0]; ii++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf]= field_line[ ii + line_size*0];
            BaryonField[Vel1Num][index_bf]= field_line[ ii + line_size*1]/field_line[ ii + line_size*0];
            BaryonField[Vel2Num][index_bf]= field_line[ ii + line_size*2]/field_line[ ii + line_size*0];
            BaryonField[Vel3Num][index_bf]= field_line[ ii + line_size*3]/field_line[ ii + line_size*0];
            CenteredB[0][index_bf]        = field_line[ ii + line_size*4];
            CenteredB[1][index_bf]        = field_line[ ii + line_size*5];
            CenteredB[2][index_bf]        = field_line[ ii + line_size*6];
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ ii + line_size*7];
              pressure[index_bf] = field_line[ii + line_size*8 ] *POW(field_line[ii + line_size*0], (Gamma - 1));
            }
           //DO collors from lines
           //DO energy from lines
         }
                //DO subgrid fluxes
                //
                //DO Magnetic fluxes (loop positions?)
          }//x sweep jj loop
      }//x sweep kk loop
    }//x sweep conditional

    //DO strang logic in here
    // Y SWEEP
    if ( 1 ) {
      for( kk = 0; kk< GridDimension[2]; kk++){
        for(ii = 0; ii<GridDimension[0]; ii++ ){
          for( jj = 0; jj<GridDimension[1]; jj++){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ jj + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*1] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*2] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*3] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*4] = CenteredB[1][index_bf];
            field_line[ jj + line_size*5] = CenteredB[2][index_bf];
            field_line[ jj + line_size*6] = CenteredB[0][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ jj + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ jj + line_size*8] = entropy[index_bf];
            }
   
            //DO fill colour lines
            //DO fill gravity line
            gravity_line[ jj ] = 0.0;
            //DO zero flux lines
          }
            //DO diffusion term
           for( jj=1; jj<GridDimension[1]; jj++){
             diffusion_line[jj] = 0.0;
           }
            FORTRAN_NAME(pde1dsolver_mhd)(field_line, &line_size, &nu, 
                &startindex, &endindex,
                CellWidthTemp[1],  &dtFixed, 
                flux_magnetic_line, //kill this
                flux_line, 
                flux_def_line, //kill this too
                diffusion_line,
                &Gamma, &csmin, &rhomin,
                &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
                &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
                &tdum0, &boxl0, &hubb, &zr, 
                &CycleNumber, &GravityOn, gravity_line, 
                a, &EquationOfState, &IsothermalSoundSpeed, &hack);


          for(jj = 0; jj<GridDimension[1]; jj++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf]= field_line[ jj + line_size*0];
            BaryonField[Vel1Num][index_bf]= field_line[ jj + line_size*3]/field_line[ jj + line_size*0];
            BaryonField[Vel2Num][index_bf]= field_line[ jj + line_size*1]/field_line[ jj + line_size*0];
            BaryonField[Vel3Num][index_bf]= field_line[ jj + line_size*2]/field_line[ jj + line_size*0];
            CenteredB[0][index_bf]        = field_line[ jj + line_size*6];
            CenteredB[1][index_bf]        = field_line[ jj + line_size*4];
            CenteredB[2][index_bf]        = field_line[ jj + line_size*5];
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ jj + line_size*7];
              pressure[index_bf] = field_line[jj + line_size*8 ] *POW(field_line[jj + line_size*0], (Gamma - 1));
            }
           //DO collors from lines
           //DO energy from lines
         }
                //DO subgrid fluxes
                //
                //DO Magnetic fluxes (loop positions?)
          }//y sweep jj loop
      }//y sweep kk loop
    }//strang conditional
    
    // z sweep
    if( 1 ){
      for(ii = 0; ii<GridDimension[0]; ii++ ){ //DO is this the fastest order?
        for( jj = 0; jj<GridDimension[1]; jj++){
          for( kk = 0; kk< GridDimension[2]; kk++){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ kk + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*1] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*2] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*3] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*4] = CenteredB[2][index_bf];
            field_line[ kk + line_size*5] = CenteredB[0][index_bf];
            field_line[ kk + line_size*6] = CenteredB[1][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ kk + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ kk + line_size*8] = entropy[index_bf];
            }
   
            //DO fill colour lines
            //DO fill gravity line
            gravity_line[ kk ] = 0.0;
            //DO zero flux lines
          }
            //DO diffusion term
           for( kk=1; kk<GridDimension[1]; kk++){
             diffusion_line[kk] = 0.0;
           }
            FORTRAN_NAME(pde1dsolver_mhd)(field_line, &line_size, &nu, 
                &startindex, &endindex,
                CellWidthTemp[2],  &dtFixed, 
                flux_magnetic_line, //kill this
                flux_line, 
                flux_def_line, //kill this too
                diffusion_line,
                &Gamma, &csmin, &rhomin,
                &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
                &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
                &tdum0, &boxl0, &hubb, &zr, 
                &CycleNumber, &GravityOn, gravity_line, 
                a, &EquationOfState, &IsothermalSoundSpeed, &hack);


          for(kk = 0; kk<GridDimension[1]; kk++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf]= field_line[ kk + line_size*0];
            BaryonField[Vel1Num][index_bf]= field_line[ kk + line_size*2]/field_line[ kk + line_size*0];
            BaryonField[Vel2Num][index_bf]= field_line[ kk + line_size*3]/field_line[ kk + line_size*0];
            BaryonField[Vel3Num][index_bf]= field_line[ kk + line_size*1]/field_line[ kk + line_size*0];
            CenteredB[0][index_bf]        = field_line[ kk + line_size*5];
            CenteredB[1][index_bf]        = field_line[ kk + line_size*6];
            CenteredB[2][index_bf]        = field_line[ kk + line_size*4];
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ kk + line_size*7];
              pressure[index_bf] = field_line[kk + line_size*8 ] *POW(field_line[kk + line_size*0], (Gamma - 1));
            }
           //DO collors from lines
           //DO energy from lines
         }
                //DO subgrid fluxes
                //
                //DO Magnetic fluxes (loop positions?)
          }//y sweep jj loop
      }//y sweep kk loop
    }

  }//strang order loop

  // DO Y LOOP
  // DO Z LOOP
  //
  //
  if( DualEnergyFormalism ){
    for( ii=0; ii<size; ii++){
      BaryonField[GENum][ii] /= BaryonField[DensNum][ii]*(Gamma-1);
    }
  }

  //DO magnetic cosmology
  //DO entropy to pressure, if applicable
  //DO clean up arrays.  Check all news in code.
  //DO not filling entropy array fails, even with DEF off.  Sort out why.
  
  return SUCCESS;
}
