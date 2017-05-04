/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR HydroSolver)
/
/  written by: Xinyu Li
/  date:       November, 2016
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     end     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     gridvel - bulk grid velocity (vector of length 3)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (not currently used)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     bottom  - true (1) if this is the lowest level
c     minsupecoef - coefficient for minimum pressure support (0 - not used)
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

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
#include "fortran.def"
#define IDX(a,b,c) ( ((c)*jn + (b))*in + (a) )


int grid::ComputeQuantumAcceleration(float *d, float dx[], float dy[], float dz[], float lapcoef)
{

  /*  Locals */

  int dim,i, j, k;

  int rank = GridRank;
  int in = GridDimension[0];
  int jn = GridDimension[1];
  int kn = GridDimension[2];

  /* Allocate temporary space */

  int size = in*jn*kn;
  float *p = new float[size];
  float *logd = new float[size];

  /* Compute log of density */
   for (i = 0; i < size; i++) {
    if (d[i] < 0) {
      ENZO_FAIL("Compute Log Density: Negative Density! \n");
    } else{
      logd[i] = log(d[i]);
      p[i] = 0;
    }
  }

    /* Loop over dimensions and difference acceleration. */
 //if( SelfGravity == 0 ){
  for (dim = 0; dim < rank; dim++) {
 
    /* Allocate acceleration field if not allocated. */
 
    if (AccelerationField[dim] == NULL)
    for (dim = 0; dim < GridRank; dim++) {
      AccelerationField[dim] = new float[size];
      for (i = 0; i < size; i++)
  AccelerationField[dim][i] = 0;
    }
  }

  if(SelfGravity==0){
    for (dim = 0; dim < rank; dim++) {

    for (i = 0; i < size; i++)
  AccelerationField[dim][i] = 0;
}
}

    /* compute Quantum Pressure Field */

    for (k = 0; k < kn; k++) {
        for (j = 0; j < jn; j++) {
            for (i = 0; i < in; i++){
          //4th order differentiation
          p[IDX(i,j,k)] = (-1./12.*(logd[IDX(max(i-2,0),j,k)]+logd[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(logd[IDX(max(i-1,0),j,k)]+logd[IDX(min(i+1,in-1),j,k)])
                    -5./2.*logd[IDX(i,j,k)])/(dx[i]*dx[i])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
                  + pow((1./12.*(logd[IDX(max(i-2,0),j,k)]-logd[IDX(min(i+2,in-1),j,k)])
                    -2./3.*(logd[IDX(max(i-1,0),j,k)]-logd[IDX(min(i+1,in-1),j,k)]))/dx[i],2)/4.;
          if (rank > 1) {
     
          p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,max(j-2,0),k)]+logd[IDX(i,min(j+2,jn-1),k)])
          				  +4./3.*(logd[IDX(i,max(j-1,0),k)]+logd[IDX(i,min(j+1,jn-1),k)])
          				  -5./2.*logd[IDX(i,j,k)])/(dy[j]*dy[j])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,max(j-2,0),k)]-logd[IDX(i,min(j+2,jn-1),k)])
          				  -2./3.*(logd[IDX(i,max(j-1,0),k)]-logd[IDX(i,min(j+1,jn-1),k)]))/dy[j],2)/4.;

          }//endif rank >1
          if (rank > 2) {
    
         p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,j,max(k-2,0))]+logd[IDX(i,j,min(k+2,kn-1))])
          				  +4./3.*(logd[IDX(i,j,max(k-1,0))]+logd[IDX(i,j,min(k+1,kn-1))])
          				  -5./2.*logd[IDX(i,j,k)])/(dz[j]*dz[j])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,j,max(k-2,0))]-logd[IDX(i,j,min(k+2,kn-1))])
          				  -2./3.*(logd[IDX(i,j,max(k-1,0))]-logd[IDX(i,j,min(k+1,kn-1))]))/dz[j],2)/4.;

          }//endif rank>2

          p[IDX(i,j,k)] = p[IDX(i,j,k)]*lapcoef;
            }//end loop over i
          }//end loop over j
        }//end loop over k

      /* Update Acceleration Field */
      /* edge-centered acceleration for Zeus, and body centered for everything else*/

      for (k = 0; k < kn; k++) {
        for (j = 0; j < jn; j++) {
          for (i = 0; i < in; i++){

            if (HydroMethod == Zeus_Hydro){
              AccelerationField[0][IDX(i,j,k)] += (p[IDX(i,j,k)]-p[IDX(max(0,i-1),j,k)])/dx[i];

            }else{
              AccelerationField[0][IDX(i,j,k)] += (p[IDX(min(in-1,i+1),j,k)]-p[IDX(max(0,i-1),j,k)])/(2*dx[i]);

            }

           

          if (rank > 1) {
     
          if (HydroMethod == Zeus_Hydro){
              AccelerationField[1][IDX(i,j,k)] += (p[IDX(i,j,k)]-p[IDX(i,max(0,j-1),k)])/dy[j];
            }else{
              AccelerationField[1][IDX(i,j,k)] += (p[IDX(i,min(jn-1,j+1),k)]-p[IDX(i,max(0,j-1),k)])/(2*dy[j]);
            }
            }//endif rank>1

          if (rank > 2) {
    
          if (HydroMethod == Zeus_Hydro){
              AccelerationField[2][IDX(i,j,k)] += (p[IDX(i,j,k)]-p[IDX(i,j,max(0,k-1))])/dz[k];
            }else{
              AccelerationField[2][IDX(i,j,k)] += (p[IDX(i,j,min(kn-1,k+1))]-p[IDX(i,j,max(0,k-1))])/(2*dz[k]);
            }

          }//endif rank>2

            }//end loop over i
          }//end loop over j
        }//end loop over k
 

  delete [] p;
  delete [] logd;
  
  
  return SUCCESS;

}
