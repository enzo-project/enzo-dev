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

int rk4(float *repsi, float *impsi,
         int in, int jn, int kn, int rank,
         float dt, float dx[], float dy[], float dz[],
         double hmcoef)

{
  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Local declarations */

  int i, j, k, jsm1, ksm1, jep1, kep1, ism2, jsm2, ksm2, jep2, kep2;
  int size = in*jn*kn;

  float *r1 = new float[size];
  float *i1 = new float[size];
  float *r2 = new float[size];
  float *i2 = new float[size];
  float *r3 = new float[size];
  float *i3 = new float[size];
  float *r4 = new float[size];
  float *i4 = new float[size];

  /* ======================================================================= */

  /* Compute varients on start indexes */

  /*jsm1 = max(js-1, 0);
  jsm2 = max(js-2, 0);
  jep1 = min(je+1, jn-1);
  jep2 = min(je+2, jn-1);
  ksm1 = max(ks-1, 0);
  ksm2 = max(ks-2, 0);
  kep1 = min(ke+1, kn-1);
  kep2 = min(ke+2, kn-1);*/

  /* Compute the laplacian term */
 // fprintf(stderr, " ksm1,kep1 %"ISYM", %"ISYM"\n", ks, ke);

  // 1st step

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          r1[IDX(i,j,k)] = -1.*(-1./12.*(impsi[IDX(max(i-2,0),j,k)]+impsi[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(impsi[IDX(max(i-1,0),j,k)]+impsi[IDX(min(i+1,in-1),j,k)])
                    -5./2.*impsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef ;
          i1[IDX(i,j,k)] = (-1./12.*(repsi[IDX(max(i-2,0),j,k)]+repsi[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(repsi[IDX(max(i-1,0),j,k)]+repsi[IDX(min(i+1,in-1),j,k)])
                    -5./2.*repsi[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef ;

          if (rank > 1) {
          r1[IDX(i,j,k)] += -1.*(-1./12.*(impsi[IDX(i,max(j-2,0),k)]+impsi[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(impsi[IDX(i,max(j-1,0),k)]+impsi[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*impsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef ;
          i1[IDX(i,j,k)] += (-1./12.*(repsi[IDX(i,max(j-2,0),k)]+repsi[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(repsi[IDX(i,max(j-1,0),k)]+repsi[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*repsi[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef ;

          }// end rank > 1
          if (rank > 2) {
          r1[IDX(i,j,k)] += -1.*(-1./12.*(impsi[IDX(i,j,max(k-2,0))]+impsi[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(impsi[IDX(i,j,max(k-1,0))]+impsi[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*impsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef ;
          i1[IDX(i,j,k)] += (-1./12.*(repsi[IDX(i,j,max(k-2,0))]+repsi[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(repsi[IDX(i,j,max(k-1,0))]+repsi[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*repsi[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef ;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */

  // 2nd step

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          r2[IDX(i,j,k)] = r1[IDX(i,j,k)] - (-1./12.*(i1[IDX(max(i-2,0),j,k)]+i1[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(i1[IDX(max(i-1,0),j,k)]+i1[IDX(min(i+1,in-1),j,k)])
                    -5./2.*i1[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt/2.;
          i2[IDX(i,j,k)] = i1[IDX(i,j,k)] + (-1./12.*(r1[IDX(max(i-2,0),j,k)]+r1[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(r1[IDX(max(i-1,0),j,k)]+r1[IDX(min(i+1,in-1),j,k)])
                    -5./2.*r1[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt/2.;

          if (rank > 1) {
          r2[IDX(i,j,k)] += -(-1./12.*(i1[IDX(i,max(j-2,0),k)]+i1[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(i1[IDX(i,max(j-1,0),k)]+i1[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*i1[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt/2.;
          i2[IDX(i,j,k)] +=  (-1./12.*(r1[IDX(i,max(j-2,0),k)]+r1[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(r1[IDX(i,max(j-1,0),k)]+r1[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*r1[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt/2.;

          }// end rank > 1
          if (rank > 2) {
          r2[IDX(i,j,k)] += -(-1./12.*(i1[IDX(i,j,max(k-2,0))]+i1[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(i1[IDX(i,j,max(k-1,0))]+i1[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*i1[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt/2.;
          i2[IDX(i,j,k)] +=  (-1./12.*(r1[IDX(i,j,max(k-2,0))]+r1[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(r1[IDX(i,j,max(k-1,0))]+r1[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*r1[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt/2.;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */

  // 3nd step

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          r3[IDX(i,j,k)] = r1[IDX(i,j,k)] - (-1./12.*(i2[IDX(max(i-2,0),j,k)]+i2[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(i2[IDX(max(i-1,0),j,k)]+i2[IDX(min(i+1,in-1),j,k)])
                    -5./2.*i2[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt/2.;
          i3[IDX(i,j,k)] = i1[IDX(i,j,k)] + (-1./12.*(r2[IDX(max(i-2,0),j,k)]+r2[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(r2[IDX(max(i-1,0),j,k)]+r2[IDX(min(i+1,in-1),j,k)])
                    -5./2.*r2[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt/2.;

          if (rank > 1) {
          r3[IDX(i,j,k)] += -(-1./12.*(i2[IDX(i,max(j-2,0),k)]+i2[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(i2[IDX(i,max(j-1,0),k)]+i2[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*i2[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt/2.;
          i3[IDX(i,j,k)] +=  (-1./12.*(r2[IDX(i,max(j-2,0),k)]+r2[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(r2[IDX(i,max(j-1,0),k)]+r2[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*r2[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt/2.;

          }// end rank > 1
          if (rank > 2) {
          r3[IDX(i,j,k)] += -(-1./12.*(i2[IDX(i,j,max(k-2,0))]+i2[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(i2[IDX(i,j,max(k-1,0))]+i2[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*i2[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt/2.;
          i3[IDX(i,j,k)] +=  (-1./12.*(r2[IDX(i,j,max(k-2,0))]+r2[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(r2[IDX(i,j,max(k-1,0))]+r2[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*r2[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt/2.;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */

  // 4th step

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          r4[IDX(i,j,k)] = r1[IDX(i,j,k)] - (-1./12.*(i3[IDX(max(i-2,0),j,k)]+i3[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(i3[IDX(max(i-1,0),j,k)]+i3[IDX(min(i+1,in-1),j,k)])
                    -5./2.*i3[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt;
          i4[IDX(i,j,k)] = i1[IDX(i,j,k)] + (-1./12.*(r3[IDX(max(i-2,0),j,k)]+r3[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(r3[IDX(max(i-1,0),j,k)]+r3[IDX(min(i+1,in-1),j,k)])
                    -5./2.*r3[IDX(i,j,k)])/(dx[i]*dx[i])/2.*hmcoef * dt;

          if (rank > 1) {
          r4[IDX(i,j,k)] += -(-1./12.*(i3[IDX(i,max(j-2,0),k)]+i3[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(i3[IDX(i,max(j-1,0),k)]+i3[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*i3[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt;
          i4[IDX(i,j,k)] +=  (-1./12.*(r3[IDX(i,max(j-2,0),k)]+r3[IDX(i,min(j+2,jn-1),k)])
                    +4./3.*(r3[IDX(i,max(j-1,0),k)]+r3[IDX(i,min(j+1,jn-1),k)])
                    -5./2.*r3[IDX(i,j,k)])/(dy[j]*dy[j])/2.*hmcoef * dt;

          }// end rank > 1
          if (rank > 2) {
          r4[IDX(i,j,k)] += -(-1./12.*(i3[IDX(i,j,max(k-2,0))]+i3[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(i3[IDX(i,j,max(k-1,0))]+i3[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*i3[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt;
          i4[IDX(i,j,k)] +=  (-1./12.*(r3[IDX(i,j,max(k-2,0))]+r3[IDX(i,j,min(k+2,kn-1))])
                    +4./3.*(r3[IDX(i,j,max(k-1,0))]+r3[IDX(i,j,min(k+1,kn-1))])
                    -5./2.*r3[IDX(i,j,k)])/(dz[k]*dz[k])/2.*hmcoef * dt;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */



  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          repsi[IDX(i,j,k)] = repsi[IDX(i,j,k)] + (r1[IDX(i,j,k)]+2.*r2[IDX(i,j,k)]+2.*r3[IDX(i,j,k)]+r4[IDX(i,j,k)])/6. *dt;
          impsi[IDX(i,j,k)] = impsi[IDX(i,j,k)] + (i1[IDX(i,j,k)]+2.*i2[IDX(i,j,k)]+2.*i3[IDX(i,j,k)]+i4[IDX(i,j,k)])/6. *dt;
      }// end loop over i
    } // end loop over j
  } // end loop over k */



delete [] r1;
delete [] i1;
delete [] r2;
delete [] i2;
delete [] r3;
delete [] i3;
delete [] r4;
delete [] i4;


  return SUCCESS;
}
