/***********************************************************************
c
c  one step evolution of RK3-TVD
c
c  written by: Xinyu Li
c  date:       May, 2017
c  modified1: 
c
c  PURPOSE:
c     Unitary operator of adding potential to Schrodinger Equation
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
#define GIDX(a,b,c) ( ((c)*(jn+2) + (b))*(in+2) + (a) )


int SchrdingerAddPotential(double *repsi, double *impsi,
         int in, int jn, int kn, int rank,
         int is, int ie, int js, int je, int ks, int ke, 
         double dt, 
         double hmcoef,
         double *p, int start1, int start2, int start3)
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

  double auxre, auxim;  

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
        auxre = cos(p[GIDX(i+start1, j+start2, k+start3)] / hmcoef * dt /2. ) * repsi[IDX(i,j,k)] 
              + sin(p[GIDX(i+start1, j+start2, k+start3)] / hmcoef * dt /2. ) * impsi[IDX(i,j,k)];

        auxim = cos(p[GIDX(i+start1, j+start2, k+start3)] / hmcoef * dt /2. ) * impsi[IDX(i,j,k)] 
              - sin(p[GIDX(i+start1, j+start2, k+start3)] / hmcoef * dt /2. ) * repsi[IDX(i,j,k)];

        repsi[IDX(i,j,k)] = auxre;
        impsi[IDX(i,j,k)] = auxim;
        //printf("potential %d %d %d\n", i+start1, j+start2, k+start3);
          

      }// end loop over i
    } // end loop over j
  } // end loop over k */

  return SUCCESS;
}
