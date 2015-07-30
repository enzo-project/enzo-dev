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

extern "C" void FORTRAN_NAME(pde1dsolver_mhd_new)(float * wx, float* colours, int * idim,int * nu,
           int * startindex, int * endindex, int * numberofcolours,
           float * dx, float * dtstrang,
           float * fluxBx, //remove
           float * fluxx, 
           float * fluxEx, //remove
           float * fluxcolour,
           float * diffcoefx, 
           float * gamma, float * csmin, float * rhomin, 
           int * MHDCTDualEnergyMethod, int * MHDCTSlopeLimiter,int * RiemannSolver, 
           int * ReconstructionMethod, int * idiffusion, int * MHDCTPowellSource,
           float * Theta,
           int * nhy, int * gravityon, float * gravityx, 
           FLOAT * a, int * EquationOfState, float * SoundSpeed, int * hack2);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

void TransverseMagneticFlux(float *Vx, float * Vy , float * Vz,
                           float *Bx, float * By , float * Bz,
                           float * F1, float * F2, int * CellDims, int* FaceDims)
{
  //For GridRank<3, the electric field still needs to be computed from transverse velocity and magnetic field.
  //This routine fills only the transverse fluxes necessary for electric field computation.
  //This is often refered to as "2.5 dimensions" in the literature.
  int index_cell, index_face, i,j,k;
  for( k=0;k<CellDims[2];k++)
    for(j=0;j<CellDims[1];j++)
      for(i=0;i<CellDims[0];i++){
        index_face = i + FaceDims[0]*(j+FaceDims[1]*k);
        index_cell = i + CellDims[0]*(j+CellDims[1]*k);
        F1[index_face] = -(Vy[index_cell]*Bx[index_cell]-Vx[index_cell]*By[index_cell]);
        F2[index_face] = -(Vz[index_cell]*Bx[index_cell]-Vx[index_cell]*Bz[index_cell]);
      }
}

int grid::SolveMHD_Li(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[],
          float ** Fluxes)
{
    

  //Big questions: 
  //1.) Strang?  Do I want to do the full pairwise strang, or the simpler style?
  //2.) What is needed for rank < 3?
  //DO clean up arrays.  Check all news in code.
  //DO check UseSpecificEnergy
  //DO strang
  //DO 2d
  
  FLOAT a[4];
  
  if(ComovingCoordinates==1){
      CosmologyComputeExpansionFactor(Time, &a[0], &a[1]) ;
      CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed, &a[2], &a[3]);
  }else{
      a[0] = 1.0;
      a[1] = 0.0;
      a[2] = 1.0;
      a[3] = 0.0;
  }


  int line_size = max( GridDimension[0], max(GridDimension[1], GridDimension[2]));
  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int ixyz = CycleNumber % 3;
  int nxz, nyz, nzz, nColour;
  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;
  int n, ii, jj, kk, index_bf,  dim, subgrid;
  int hack = 0; //a flag for testing the solver.
  float dtdx;

  int nu = 6;
  float csmin = 1e-13, rhomin = 1e-6;

  int line_width = 9;  //the number of conserved quantities.
  float * field_line     = new float[line_size * line_width];
  float * flux_line      = new float[line_size * line_width];
  float * colour_line    = new float[line_size * NumberOfColours];
  float * flux_colour    = new float[line_size * NumberOfColours];
  float * gravity_line   = new float[line_size];
  float * diffusion_line = new float[line_size];

  //Both of these need to be convolved with the above.
  float * flux_magnetic_line = new float[ line_size ]; //for Powell fluxes
  float * flux_def_line = new float[ line_size ]; //For dual energy.  Get rid of this.
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  float inv_sqrt_a = 1./sqrt(a[0]), sqrt_a = sqrt(a[0]);
  if( ComovingCoordinates ){
    if ( EquationOfState == 0 ){
      for( ii=0; ii<size; ii++){
        BaryonField[TENum][ii] -= 0.5*(BaryonField[B1Num][ii]*BaryonField[B1Num][ii]+
                                       BaryonField[B2Num][ii]*BaryonField[B2Num][ii]+
                                       BaryonField[B3Num][ii]*BaryonField[B3Num][ii])/BaryonField[DensNum][ii];
      }
    }

    for( ii=0; ii<size; ii++){
      BaryonField[B1Num][ii] *= inv_sqrt_a;
      BaryonField[B2Num][ii] *= inv_sqrt_a;
      BaryonField[B3Num][ii] *= inv_sqrt_a;
    }
    
    if ( EquationOfState == 0 ){
      for( ii=0; ii<size; ii++){
        BaryonField[TENum][ii] += 0.5*(BaryonField[B1Num][ii]*BaryonField[B1Num][ii]+
                                       BaryonField[B2Num][ii]*BaryonField[B2Num][ii]+
                                       BaryonField[B3Num][ii]*BaryonField[B3Num][ii])/BaryonField[DensNum][ii];
      }
    }
  }//comoving

 

  float * pressure;
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
            -0.5*(BaryonField[B1Num][ii]*BaryonField[B1Num][ii]+
                  BaryonField[B2Num][ii]*BaryonField[B2Num][ii]+
                  BaryonField[B3Num][ii]*BaryonField[B3Num][ii])); 
      }
    }
  }

  //Pointers to magnetic fluxes for simplicity below.  
  float * MagFluxX1 = Fluxes[0],
        * MagFluxX2 = Fluxes[0]+MagneticSize[0],
        * MagFluxY1 = Fluxes[1],
        * MagFluxY2 = Fluxes[1]+MagneticSize[1],
        * MagFluxZ1 = Fluxes[2],
        * MagFluxZ2 = Fluxes[2]+MagneticSize[2];

  int *fistart[3], *fiend[3], *fjstart[3], *fjend[3], *nfi[3], *lindex[3], *rindex[3];
  for ( dim=0;dim<3;dim++){
    fistart[dim] = new int[NumberOfSubgrids];
    fiend[dim]   = new int[NumberOfSubgrids];
    fjstart[dim] = new int[NumberOfSubgrids];
    fjend[dim]   = new int[NumberOfSubgrids];
    nfi[dim]     = new int[NumberOfSubgrids];
    lindex[dim]  = new int[NumberOfSubgrids];
    rindex[dim]  = new int[NumberOfSubgrids];

  }

  int idim, jdim, offset;
  for( n=0; n<NumberOfSubgrids; n++){
    //Transverse start and end index (y and z coordinates for x flux, etc.)
    for( dim=0;dim<3;dim++){
      idim = (dim == 0 ) ? 1 : ( ( dim == 1 ) ? 0 : 0 );
      jdim = (dim == 0 ) ? 2 : ( ( dim == 1 ) ? 2 : 1 );
      fistart[dim][n] = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] - GridGlobalStart[idim];
      fiend[dim][n]   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim]   - GridGlobalStart[idim];
      fjstart[dim][n] = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - GridGlobalStart[jdim];
      fjend[dim][n]   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim]   - GridGlobalStart[jdim];
      //width of plane for computing offsets
      nfi[dim][n]     = fiend[dim][n] - fistart[dim][n] + 1;
      lindex[dim][n] = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] - GridGlobalStart[dim]-1; //-1 for fortran indices
      rindex[dim][n] = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] - GridGlobalStart[dim];   //-1 for fortran, +1 for the right face.
    }
    //Position of planes, longitudinal (x position for x flux, etc)
  }//index allocation

  //strang loop.
  
  int startindex, endindex;
  for (n = ixyz; n < ixyz+3; n++) {
    startindex = 3;
    endindex = GridDimension[0] - 2;
    dtdx = dtFixed/CellWidthTemp[0][0];
    //x loope
    if ( (n % 3 == 0) && nxz > 1){ 
      for( kk = 0; kk< GridDimension[2]; kk++){
        for( jj = 0; jj<GridDimension[1]; jj++){
          for(ii = 0; ii<GridDimension[0]; ii++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ ii + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*1] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*2] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*3] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ ii + line_size*4] = BaryonField[B1Num][index_bf];
            field_line[ ii + line_size*5] = BaryonField[B2Num][index_bf];
            field_line[ ii + line_size*6] = BaryonField[B3Num][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ ii + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ ii + line_size*8] = pressure[index_bf]/POW(BaryonField[DensNum][index_bf],Gamma-1);
            }
            if( GravityOn )
              gravity_line[ ii ] = AccelerationField[0][index_bf];
            if( NumberOfColours > 0){
              for( nColour=0; nColour<NumberOfColours; nColour++){
                colour_line[ ii + line_size*nColour ] = BaryonField[ colnum[nColour] ][index_bf];
              }
            }//colours?
          }//ii loop
          //DO diffusion term
          for( ii=1; ii<GridDimension[0]; ii++){
            diffusion_line[ii] = 0.0;
          }
          FORTRAN_NAME(pde1dsolver_mhd_new)(field_line, colour_line, &line_size, &nu, 
            &startindex, &endindex, &NumberOfColours,
            CellWidthTemp[0],  &dtFixed, 
            flux_magnetic_line, //kill this
            flux_line, 
            flux_def_line, //kill this too
            flux_colour,
            diffusion_line,
            &Gamma, &csmin, &rhomin,
            &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
            &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
            &Theta_Limiter,
            &CycleNumber, &GravityOn, gravity_line, 
            a, &EquationOfState, &IsothermalSoundSpeed, &hack);

          for(ii = 0; ii<GridDimension[0]; ii++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf] = field_line[ ii + line_size*0];
            BaryonField[Vel1Num][index_bf] = field_line[ ii + line_size*1]/field_line[ ii + line_size*0];
            BaryonField[Vel2Num][index_bf] = field_line[ ii + line_size*2]/field_line[ ii + line_size*0];
            BaryonField[Vel3Num][index_bf] = field_line[ ii + line_size*3]/field_line[ ii + line_size*0];
            BaryonField[B1Num][index_bf]         = field_line[ ii + line_size*4];
            BaryonField[B2Num][index_bf]         = field_line[ ii + line_size*5];
            BaryonField[B3Num][index_bf]         = field_line[ ii + line_size*6];
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ ii + line_size*7];
              pressure[index_bf]           = field_line[ii + line_size*8 ] *POW(field_line[ii + line_size*0], (Gamma - 1));
            }
            if( NumberOfColours > 0 ){
              for( nColour=0; nColour<NumberOfColours; nColour++){
                BaryonField[colnum[nColour]][index_bf] = colour_line[ ii + line_size*nColour];
              }
            }//colours

          }//ii
          //Fill subgrids.
          dim = 0;
          //xflux
          for( subgrid=0; subgrid<NumberOfSubgrids; subgrid++){
            if( ( jj >= fistart[dim][subgrid] && jj <= fiend[dim][subgrid] ) &&
                ( kk >= fjstart[dim][subgrid] && kk <= fjend[dim][subgrid] ) ){
              offset = (jj - fistart[dim][subgrid] ) + (kk - fjstart[dim][subgrid] )*nfi[dim][subgrid];
              ii = lindex[dim][subgrid];
              SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim][offset] = dtdx*flux_line[ii + line_size*0];
              SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim][offset] = dtdx*flux_line[ii + line_size*1];
              SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim][offset] = dtdx*flux_line[ii + line_size*2];
              SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim][offset] = dtdx*flux_line[ii + line_size*3];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->LeftFluxes[TENum][dim][offset] = dtdx*flux_line[ii + line_size*7];
              }
              if( DualEnergyFormalism ){
                SubgridFluxes[subgrid]->LeftFluxes[GENum][dim][offset] = dtdx*flux_def_line[ii];
              }
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->LeftFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[ii +line_size*nColour];
                }
              }
              ii = rindex[dim][subgrid];
              SubgridFluxes[subgrid]->RightFluxes[DensNum][dim][offset]= dtdx*flux_line[ii + line_size*0];
              SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim][offset]= dtdx*flux_line[ii + line_size*1];
              SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim][offset]= dtdx*flux_line[ii + line_size*2];
              SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim][offset]= dtdx*flux_line[ii + line_size*3];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->RightFluxes[TENum][dim][offset]= dtdx*flux_line[ii + line_size*7];
              }
              if( DualEnergyFormalism ){
                SubgridFluxes[subgrid]->RightFluxes[GENum][dim][offset]= dtdx*flux_def_line[ii];
              }//GE flux
              //color fluxes
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->RightFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[ii +line_size*nColour];
                }
              }
                //end color fluxes
              
            }//subgrid ok.
          }//subgrid loop
          for(ii = 3; ii<GridDimension[0]-1; ii++ ){
            index_bf = ii + MagneticDims[0][0]*(jj + MagneticDims[0][1]*kk);
            MagFluxX1[index_bf] = flux_line[ ii-1 + line_size*5];
            MagFluxX2[index_bf] = flux_line[ ii-1 + line_size*6];
          }
        }//x sweep jj loop
      }//x sweep kk loop
    }//x sweep conditional
    
    startindex = 3;
    endindex = GridDimension[1] - 2;
    dtdx = dtFixed/CellWidthTemp[1][0];
    // Y SWEEP
    if ( (n % 3 == 1) )
      if ( nyz == 1 ){
        TransverseMagneticFlux(BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[Vel1Num],
                               BaryonField[B2Num], BaryonField[B3Num], BaryonField[B1Num], //CenteredB[1], CenteredB[2], CenteredB[0], 
                               MagFluxY2, MagFluxY1, GridDimension, MagneticDims[1]);

      }else if ( nyz > 1) {
      for( kk = 0; kk< GridDimension[2]; kk++){
        for(ii = 0; ii<GridDimension[0]; ii++ ){
          for( jj = 0; jj<GridDimension[1]; jj++){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ jj + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*1] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*2] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*3] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ jj + line_size*4] = BaryonField[B2Num][index_bf];// CenteredB[1][index_bf];
            field_line[ jj + line_size*5] = BaryonField[B3Num][index_bf];// CenteredB[2][index_bf];
            field_line[ jj + line_size*6] = BaryonField[B1Num][index_bf];// CenteredB[0][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ jj + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ jj + line_size*8] = pressure[index_bf]/POW(BaryonField[DensNum][index_bf],Gamma-1);
            }
            for( nColour=0; nColour<NumberOfColours; nColour++){
              colour_line[ jj + line_size*nColour ] = BaryonField[ colnum[nColour] ][index_bf];
            }
   
            if( GravityOn )
              gravity_line[ jj ] = AccelerationField[1][index_bf];
          }
          for( jj=1; jj<GridDimension[1]; jj++){
             diffusion_line[jj] = 0.0;
          }
          FORTRAN_NAME(pde1dsolver_mhd_new)(field_line, colour_line, &line_size, &nu, 
            &startindex, &endindex, &NumberOfColours,
            CellWidthTemp[1],  &dtFixed, 
            flux_magnetic_line, //kill this
            flux_line, 
            flux_def_line, //kill this too
            flux_colour,
            diffusion_line,
            &Gamma, &csmin, &rhomin,
            &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
            &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
            &Theta_Limiter,
            &CycleNumber, &GravityOn, gravity_line, 
            a, &EquationOfState, &IsothermalSoundSpeed, &hack);


          for(jj = 0; jj<GridDimension[1]; jj++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf]= field_line[ jj + line_size*0];
            BaryonField[Vel1Num][index_bf]= field_line[ jj + line_size*3]/field_line[ jj + line_size*0];
            BaryonField[Vel2Num][index_bf]= field_line[ jj + line_size*1]/field_line[ jj + line_size*0];
            BaryonField[Vel3Num][index_bf]= field_line[ jj + line_size*2]/field_line[ jj + line_size*0];
            BaryonField[B1Num][index_bf]       = field_line[ jj + line_size*6];//           CenteredB[0][index_bf] 
            BaryonField[B2Num][index_bf]       = field_line[ jj + line_size*4];//           CenteredB[1][index_bf] 
            BaryonField[B3Num][index_bf]       = field_line[ jj + line_size*5];//           CenteredB[2][index_bf] 
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ jj + line_size*7];
              pressure[index_bf] = field_line[jj + line_size*8 ] *POW(field_line[jj + line_size*0], (Gamma - 1));
            }
            if( NumberOfColours > 0 ){
              for( nColour=0; nColour<NumberOfColours; nColour++){
                BaryonField[colnum[nColour]][index_bf] = colour_line[ jj + line_size*nColour];
              }
            }//colours
          }
          dim = 1;
          //yflux
          for( subgrid=0; subgrid<NumberOfSubgrids; subgrid++){
            if( ( ii >= fistart[dim][subgrid] && ii <= fiend[dim][subgrid] ) &&
                ( kk >= fjstart[dim][subgrid] && kk <= fjend[dim][subgrid] ) ){
              offset = (ii - fistart[dim][subgrid] ) + (kk - fjstart[dim][subgrid] )*nfi[dim][subgrid];
              jj = lindex[dim][subgrid];
              SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim][offset] = dtdx*flux_line[jj + line_size*0];
              SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim][offset] = dtdx*flux_line[jj + line_size*3];
              SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim][offset] = dtdx*flux_line[jj + line_size*1];
              SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim][offset] = dtdx*flux_line[jj + line_size*2];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->LeftFluxes[TENum][dim][offset]   = dtdx*flux_line[jj + line_size*7];
              }//te flux
              if( DualEnergyFormalism ){
               SubgridFluxes[subgrid]->LeftFluxes[GENum][dim][offset] = dtdx*flux_def_line[jj];
              }//ge flux
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->LeftFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[jj +line_size*nColour];
                }
              }
              jj = rindex[dim][subgrid];
              SubgridFluxes[subgrid]->RightFluxes[DensNum][dim][offset]= dtdx*flux_line[jj + line_size*0];
              SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim][offset]= dtdx*flux_line[jj + line_size*3];
              SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim][offset]= dtdx*flux_line[jj + line_size*1];
              SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim][offset]= dtdx*flux_line[jj + line_size*2];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->RightFluxes[TENum][dim][offset]  = dtdx*flux_line[jj + line_size*7];
              }//TE flux
              if( DualEnergyFormalism ){
               SubgridFluxes[subgrid]->RightFluxes[GENum][dim][offset]= dtdx*flux_def_line[jj];
              }//GE flux
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->RightFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[jj +line_size*nColour];
                }
              }
              
            }//subgrid ok.
          }//subgrid loop
          for(jj = 3; jj<GridDimension[1]-1; jj++ ){
            index_bf = ii + MagneticDims[1][0]*(jj + MagneticDims[1][1]*kk);
            MagFluxY2[index_bf] = flux_line[ jj-1 + line_size*5];
            MagFluxY1[index_bf] = flux_line[ jj-1 + line_size*6];
          }
        }//y sweep ii loop
      }//y sweep kk loop
    }//strang conditional
    
    startindex = 3;
    endindex = GridDimension[2] - 2;
    dtdx = dtFixed/CellWidthTemp[2][0];
    // z sweep
    if( (n % 3 == 2) )
      if ( nzz == 1 ){
        TransverseMagneticFlux(BaryonField[Vel3Num], BaryonField[Vel1Num], BaryonField[Vel2Num],
                               BaryonField[B3Num], BaryonField[B1Num], BaryonField[B2Num],// BaryonField[BCenteredB[2], CenteredB[0], CenteredB[1], 
                               MagFluxZ1, MagFluxZ2, GridDimension, MagneticDims[2]);

      }else if ( nzz > 1) {

      for(ii = 0; ii<GridDimension[0]; ii++ ){ //DO is this the fastest order?
        for( jj = 0; jj<GridDimension[1]; jj++){
          for( kk = 0; kk< GridDimension[2]; kk++){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            field_line[ kk + line_size*0] = BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*1] = BaryonField[Vel3Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*2] = BaryonField[Vel1Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*3] = BaryonField[Vel2Num][index_bf]*BaryonField[DensNum][index_bf];
            field_line[ kk + line_size*4] = BaryonField[B3Num][index_bf];//     CenteredB[2][index_bf];
            field_line[ kk + line_size*5] = BaryonField[B1Num][index_bf];//     CenteredB[0][index_bf];
            field_line[ kk + line_size*6] = BaryonField[B2Num][index_bf];//     CenteredB[1][index_bf];
            if ( EquationOfState == 0 ){
              field_line[ kk + line_size*7] = BaryonField[TENum][index_bf];
              field_line[ kk + line_size*8] = pressure[index_bf]/POW(BaryonField[DensNum][index_bf],Gamma-1);
            }
            if( NumberOfColours > 0){
              for( nColour=0; nColour<NumberOfColours; nColour++){
                colour_line[ kk + line_size*nColour ] = BaryonField[ colnum[nColour] ][index_bf];
              }
            }//colours?
   
            if( GravityOn )
              gravity_line[ kk ] = AccelerationField[2][index_bf];
          }
            //DO diffusion term
          for( kk=1; kk<GridDimension[1]; kk++){
            diffusion_line[kk] = 0.0;
          }
          
          FORTRAN_NAME(pde1dsolver_mhd_new)(field_line, colour_line, &line_size, &nu, 
            &startindex, &endindex, &NumberOfColours,
            CellWidthTemp[2],  &dtFixed, 
            flux_magnetic_line, //kill this
            flux_line, 
            flux_def_line, //kill this too
            flux_colour,
            diffusion_line,
            &Gamma, &csmin, &rhomin,
            &MHDCTDualEnergyMethod, &MHDCTSlopeLimiter, &RiemannSolver, 
            &ReconstructionMethod, &PPMDiffusionParameter, &MHDCTPowellSource,
            &Theta_Limiter,
            &CycleNumber, &GravityOn, gravity_line, 
            a, &EquationOfState, &IsothermalSoundSpeed, &hack);
            


          for(kk = 0; kk<GridDimension[2]; kk++ ){
            index_bf = ii + GridDimension[0]*(jj + GridDimension[1]*kk);
            BaryonField[DensNum][index_bf]= field_line[ kk + line_size*0];
            BaryonField[Vel1Num][index_bf]= field_line[ kk + line_size*2]/field_line[ kk + line_size*0];
            BaryonField[Vel2Num][index_bf]= field_line[ kk + line_size*3]/field_line[ kk + line_size*0];
            BaryonField[Vel3Num][index_bf]= field_line[ kk + line_size*1]/field_line[ kk + line_size*0];
             BaryonField[B1Num][index_bf]       = field_line[ kk + line_size*5];//CenteredB[0][index_bf] 
             BaryonField[B2Num][index_bf]       = field_line[ kk + line_size*6];//CenteredB[1][index_bf] 
             BaryonField[B3Num][index_bf]       = field_line[ kk + line_size*4];//CenteredB[2][index_bf] 
            if ( EquationOfState == 0 ){
              BaryonField[TENum][index_bf] = field_line[ kk + line_size*7];
              pressure[index_bf] = field_line[kk + line_size*8 ] *POW(field_line[kk + line_size*0], (Gamma - 1));
            }
            if( NumberOfColours > 0 ){
              for( nColour=0; nColour<NumberOfColours; nColour++){
                BaryonField[colnum[nColour]][index_bf] = colour_line[ kk + line_size*nColour];
              }
            }//colours
          }
          dim = 2;
          //zflux
          for( subgrid=0; subgrid<NumberOfSubgrids; subgrid++){
            if( ( ii >= fistart[dim][subgrid] && ii <= fiend[dim][subgrid] ) &&
                ( jj >= fjstart[dim][subgrid] && jj <= fjend[dim][subgrid] ) ){
              offset = (ii - fistart[dim][subgrid] ) + (jj - fjstart[dim][subgrid] )*nfi[dim][subgrid];
              kk = lindex[dim][subgrid];
              SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim][offset] = dtdx*flux_line[kk + line_size*0];
              SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim][offset] = dtdx*flux_line[kk + line_size*2];
              SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim][offset] = dtdx*flux_line[kk + line_size*3];
              SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim][offset] = dtdx*flux_line[kk + line_size*1];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->LeftFluxes[TENum][dim][offset] = dtdx*flux_line[kk + line_size*7];
              }//te flux
              if( DualEnergyFormalism ){
                SubgridFluxes[subgrid]->LeftFluxes[GENum][dim][offset] = dtdx*flux_def_line[kk];
              }//ge flux
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->LeftFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[kk +line_size*nColour];
                }
              }
              kk = rindex[dim][subgrid];
              SubgridFluxes[subgrid]->RightFluxes[DensNum][dim][offset]= dtdx*flux_line[kk + line_size*0];
              SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim][offset]= dtdx*flux_line[kk + line_size*2];
              SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim][offset]= dtdx*flux_line[kk + line_size*3];
              SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim][offset]= dtdx*flux_line[kk + line_size*1];
              if( EquationOfState == 0 ){
                SubgridFluxes[subgrid]->RightFluxes[TENum][dim][offset]  = dtdx*flux_line[kk + line_size*7];
              }//TE flux
              
              if( DualEnergyFormalism ){
                SubgridFluxes[subgrid]->RightFluxes[GENum][dim][offset]= dtdx*flux_def_line[kk];
              }//GE flux
              if( NumberOfColours > 0 ){
                for( nColour=0; nColour<NumberOfColours; nColour++){
                  SubgridFluxes[subgrid]->RightFluxes[colnum[nColour]][dim][offset]= dtdx*flux_colour[kk +line_size*nColour];
                }
              }
            }//subgrid ok.
          }//subgrid loop
          for(kk = 3; kk<GridDimension[2]-1; kk++ ){
            index_bf = ii + MagneticDims[2][0]*(jj + MagneticDims[2][1]*kk);
            MagFluxZ1[index_bf] = flux_line[ kk-1 + line_size*5];
            MagFluxZ2[index_bf] = flux_line[ kk-1 + line_size*6];
          }
          }//z sweep jj loop
      }//z sweep kk loop
    }

  }//strang order loop

  if( DualEnergyFormalism ){
    for( ii=0; ii<size; ii++){
      BaryonField[GENum][ii] /= BaryonField[DensNum][ii]*(Gamma-1);
    }
  }

  delete [] field_line;
  delete [] flux_line;
  delete [] colour_line;
  delete [] flux_colour;
  delete [] gravity_line;
  delete [] diffusion_line;

  //Both of these need to be convolved with the above.
  delete [] flux_magnetic_line;
  delete [] flux_def_line;

  if ( ! DualEnergyFormalism )
    delete [] pressure;

  for( dim=0;dim<3;dim++){
    delete [] fistart[dim];
    delete [] fiend[dim];
    delete [] fjstart[dim];
    delete [] fjend[dim];
    delete [] nfi[dim];
    delete [] lindex[dim];
    delete [] rindex[dim];
  }

  if( ComovingCoordinates ){
    if ( EquationOfState == 0 ){
      for( ii=0; ii<size; ii++){
        BaryonField[TENum][ii] -= 0.5*(BaryonField[B1Num][ii]*BaryonField[B1Num][ii]+
                                       BaryonField[B2Num][ii]*BaryonField[B2Num][ii]+
                                       BaryonField[B3Num][ii]*BaryonField[B3Num][ii])/BaryonField[DensNum][ii];
      }
    }

    for( int field=0; field<3; field++){
      for( ii=0; ii<size; ii++){
        BaryonField[B1Num][ii] *= sqrt_a;
        BaryonField[B2Num][ii] *= sqrt_a;
        BaryonField[B3Num][ii] *= sqrt_a;
      }
      for( ii=0; ii<MagneticSize[field]; ii++){
        Fluxes[field][ii] *= sqrt_a/a[2];
        Fluxes[field][ii+MagneticSize[field]] *= sqrt_a/a[2];
      }
    }

    if ( EquationOfState == 0 ){
      for( ii=0; ii<size; ii++){
        BaryonField[TENum][ii] += 0.5*(BaryonField[B1Num][ii]*BaryonField[B1Num][ii]+
                                       BaryonField[B2Num][ii]*BaryonField[B2Num][ii]+
                                       BaryonField[B3Num][ii]*BaryonField[B3Num][ii])/BaryonField[DensNum][ii];
      }
    }
  }//comoving


  
  return SUCCESS;
}
