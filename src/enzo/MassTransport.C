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

int MassTransport(double *d, double *repsi, double *impsi,
         int in, int jn, int kn, int rank,
         int is, int ie, int js, int je, int ks, int ke, 
         double dt, double dx[], double dy[], double dz[],
         double hmcoef)

{
  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Local declarations */

  int i, j, k, jsm1, ksm1, jep1, kep1, ism2, jsm2, ksm2, jep2, kep2;
  int size = in*jn*kn;

  double *px = new double[size];
  double *py = new double[size];
  double *pz = new double[size];

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

  // 1. Modify Psi that |Psi|^2 agree with density

  /*for (i = 0; i < size; i++){
      double NormalizeFactor = d[i]/(repsi[i]*repsi[i] + impsi[i]*impsi[i]);
      repsi[i] = repsi[i]*NormalizeFactor;
      impsi[i] = impsi[i]*NormalizeFactor;
    }*/ 

  // 2. Compute Momentum Density at the Cell edge

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          px[IDX(i,j,k)] =  (repsi[IDX(i,j,k)]+repsi[IDX(max(i-1,0),j,k)])/2. 
                          * (impsi[IDX(i,j,k)]-impsi[IDX(max(i-1,0),j,k)])/dx[i]*hmcoef
                          - (impsi[IDX(i,j,k)]+impsi[IDX(max(i-1,0),j,k)])/2. 
                          * (repsi[IDX(i,j,k)]-repsi[IDX(max(i-1,0),j,k)])/dx[i]*hmcoef;

          if (rank > 1) {
          py[IDX(i,j,k)] =  (repsi[IDX(i,j,k)]+repsi[IDX(i,max(j-1,0),k)])/2. 
                          * (impsi[IDX(i,j,k)]-impsi[IDX(i,max(j-1,0),k)])/dy[j]*hmcoef
                          - (impsi[IDX(i,j,k)]+impsi[IDX(i,max(j-1,0),k)])/2. 
                          * (repsi[IDX(i,j,k)]-repsi[IDX(i,max(j-1,0),k)])/dy[j]*hmcoef;
          }// end rank > 1
          if (rank > 2) {
          pz[IDX(i,j,k)] =  (repsi[IDX(i,j,k)]+repsi[IDX(i,j,max(k-1,0))])/2 
                          * (impsi[IDX(i,j,k)]-impsi[IDX(i,j,max(k-1,0))])/dz[k]*hmcoef
                          - (impsi[IDX(i,j,k)]+impsi[IDX(i,j,max(k-1,0))])/2 
                          * (repsi[IDX(i,j,k)]-repsi[IDX(i,j,max(k-1,0))])/dz[k]*hmcoef;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */

  // 3. Update Density

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
          d[IDX(i,j,k)] = d[IDX(i,j,k)] - (px[IDX(min(i+1,in-1),j,k)] - px[IDX(i,j,k)])/dx[i]*dt;

          if (rank > 1) {
          d[IDX(i,j,k)] = d[IDX(i,j,k)] - (py[IDX(i,min(j+1,jn-1),k)] - py[IDX(i,j,k)])/dy[j]*dt;
          }// end rank > 1
          if (rank > 2) {
          d[IDX(i,j,k)] = d[IDX(i,j,k)] - (pz[IDX(i,j,min(k+1,kn-1))] - pz[IDX(i,j,k)])/dz[k]*dt;

          }// end rank > 2
      }// end loop over i
    } // end loop over j
  } // end loop over k */

delete [] px;
delete [] py;
delete [] pz;


  return SUCCESS;
}
