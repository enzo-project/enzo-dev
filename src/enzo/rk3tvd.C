/***********************************************************************
c
c  one step evolution of RK3-TVD
c
c  written by: Xinyu Li
c  date:       May, 2017
c  modified1: 
c
c  PURPOSE:
c     Adds the source substeps
c
c  EXTERNALS:
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     i,j,kn  - dimensions of field arrays
c     igamfield - indicates if gamma should be a field
c     ipresfree - pressure-free flag (0 = off, 1 = on, i.e. p=0)
c     rank    - dimension of problem (not currently used)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     C1,C2   - Linear and quadratic artifificla viscosity parameters
c     minsupecoef - coefficient for minimum pressure support
c
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "fortran.def"

#define IDX(a,b,c) ( ((c)*jn + (b))*in + (a) )

int rk3tvd(double *repsi, double *impsi, double *repsi0, double *impsi0,
         double alpha, double beta,
         int in, int jn, int kn, int rank,
         int is, int ie, int js, int je, int ks, int ke, 
         double dt, double dx[], double dy[], double dz[],
         double hmcoef,
         int gravity, double *p)
{
  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Local declarations */

  int i, j, k, jsm1, ksm1, jep1, kep1, ism2, jsm2, ksm2, jep2, kep2;
  int size = in*jn*kn;

  float *delre = new float[size];
  float *delim = new float[size];

  /* ======================================================================= */

  /* Compute varients on start indexes */

  jsm1 = max(js-1, 0);
  jsm2 = max(js-2, 0);
  jep1 = min(je+1, jn-1);
  jep2 = min(je+2, jn-1);
  ksm1 = max(ks-1, 0);
  ksm2 = max(ks-2, 0);
  kep1 = min(ke+1, kn-1);
  kep2 = min(ke+2, kn-1);

  /* Compute the laplacian term */
 // fprintf(stderr, " ksm1,kep1 %"ISYM", %"ISYM"\n", ks, ke);

  

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          /*delre[IDX(i,j,k)] = -1.*((impsi[IDX(max(i-1,0),j,k)]+impsi[IDX(min(i+1,in-1),j,k)])
                    -2.*impsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef;
          delim[IDX(i,j,k)] = ((repsi[IDX(max(i-1,0),j,k)]+repsi[IDX(min(i+1,in-1),j,k)])
                    -2.*repsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef;*/

          //4th order differentiation
          delre[IDX(i,j,k)] = -1.*(-1./12.*(impsi[IDX(max(i-2,0),j,k)]+impsi[IDX(min(i+2,in-1),j,k)])
          				  +4./3.*(impsi[IDX(max(i-1,0),j,k)]+impsi[IDX(min(i+1,in-1),j,k)])
          				  -5./2.*impsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef;

          delim[IDX(i,j,k)] = (-1./12.*(repsi[IDX(max(i-2,0),j,k)]+repsi[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(repsi[IDX(max(i-1,0),j,k)]+repsi[IDX(min(i+1,in-1),j,k)])
                    -5./2.*repsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef;

          if (rank > 1) {
          /*delre[IDX(i,j,k)] += -1.*((impsi[IDX(i,max(j-1,0),k)]+impsi[IDX(i, min(j+1,jn-1),k)])
                    -2.*impsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef;
          delim[IDX(i,j,k)] += ((repsi[IDX(i, max(j-1,0),k)]+repsi[IDX(i, min(j+1,jn-1),k)])
                    -2.*repsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef;*/
          	//4th order differentiation
          delre[IDX(i,j,k)] += -1.*(-1./12.*(impsi[IDX(i,max(j-2,0),k)]+impsi[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(impsi[IDX(i,max(j-1,0),k)]+impsi[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*impsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef;

          delim[IDX(i,j,k)] += (-1./12.*(repsi[IDX(i,max(j-2,0),k)]+repsi[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(repsi[IDX(i,max(j-1,0),k)]+repsi[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*repsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef;

          }// end rank > 1
          if (rank > 2) {
          /*delre[IDX(i,j,k)] += -1.*((impsi[IDX(i, j, max(k-1,0))]+impsi[IDX(i, j, min(k+1,kn-1))])
                    -2.*impsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef;
          delim[IDX(i,j,k)] += ((repsi[IDX(i, j, max(k-1,0))]+repsi[IDX(i, j, min(k+1,kn-1))])
                    -2.*repsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef;*/
          //4th order differentiation
          delre[IDX(i,j,k)] += -1.*(-1./12.*(impsi[IDX(i,j,max(k-2,0))]+impsi[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(impsi[IDX(i,j,max(k-1,0))]+impsi[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*impsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef;

          delim[IDX(i,j,k)] += (-1./12.*(repsi[IDX(i,j,max(k-2,0))]+repsi[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(repsi[IDX(i,j,max(k-1,0))]+repsi[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*repsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef;

          }// end rank > 2

          /*if(gravity){
            delre[IDX(i,j,k)] += p[IDX(i,j,k)] * impsi[IDX(i,j,k)] / hmcoef;            
            delim[IDX(i,j,k)] -= p[IDX(i,j,k)] * repsi[IDX(i,j,k)] / hmcoef;
          }*/

      }// end loop over i
    } // end loop over j
  } // end loop over k */

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          repsi[IDX(i,j,k)] = alpha*repsi0[IDX(i,j,k)] + beta*(repsi[IDX(i,j,k)]+delre[IDX(i,j,k)]*dt);
          impsi[IDX(i,j,k)] = alpha*impsi0[IDX(i,j,k)] + beta*(impsi[IDX(i,j,k)]+delim[IDX(i,j,k)]*dt);
      }// end loop over i
    } // end loop over j
  } // end loop over k */



delete [] delre;
delete [] delim;

  return SUCCESS;
}
