/***********************************************************************
c
c  SOURCE TERMS FOR ZEUS HYDRO (CARTESIAN ONLY)
c
c  written by: Greg Bryan (implemented from Stone & Norman, ApJS 80, 753)
c  date:       February, 1997
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

int ZeusFDM(float *d, float *e, float *u, float *v, float *w, float *p,
	       int in, int jn, int kn, int rank,
	       int is, int ie, int js, int je, int ks, int ke, 
	       float C1, float C2, float *gamma, float dt, float dx[], float dy[], float dz[],
	       int gravity, float *gr_xacc, float *gr_yacc, float *gr_zacc, 
	       float minsupecoef, float lapcoef)
{
  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Local declarations */

  int i, j, k, jsm1, ksm1, jep1, kep1, ism2, jsm2, ksm2, jep2, kep2, n, nstep;
  int size = in*jn*kn;
  float gamma1 = gamma[0];
  float  deltav, dt1, e1;
  float *logd = new float[size];
  float q[ijk];

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

  /* compute log of density */
  
   for (i = 0; i < size; i++) {
    if (d[i] < 0) {
      fprintf(stderr, "u,v,w,d,e=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM"  dt=%"GSYM"\n", 
        u[i],v[i],w[i],d[i],p[i], dx[0], dt);
      ENZO_FAIL("ZeusSolver: Negative Density! \n");
    } else{
      logd[i] = log(d[i]);
    }
  }

  /* Compute the quantum pressure (this could be optimized) */

#define QP_4TH_ORDER

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){

	  //1st order differentiaion - keep just in case it is useful to explore this
#ifdef QP_1ST_ORDER
	p[IDX(i,j,k)] = (pow(d[IDX(max(i-1,0),j,k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(min(i+1,in-1),j,k)],0.5))/pow(dx[i],2);
#endif

          //2nd order differentiation - keep just in case it is useful to explore this
#ifdef QP_2ND_ORDER
	p[IDX(i,j,k)] = (logd[IDX(max(i-1,0),j,k)]-2.0*logd[IDX(i,j,k)]+logd[IDX(min(i+1,in-1),j,k)])/(dx[i]*dx[i])/2.;
	p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(min(i+1,in-1),j,k)]-logd[IDX(max(i-1,0),j,k)])/(dx[i]*2.),2)/4.;
#endif

          //4th order differentiation - use this
#ifdef QP_4TH_ORDER
          p[IDX(i,j,k)] = (-1./12.*(logd[IDX(max(i-2,0),j,k)]+logd[IDX(min(i+2,in-1),j,k)])
          				  +4./3.*(logd[IDX(max(i-1,0),j,k)]+logd[IDX(min(i+1,in-1),j,k)])
          				  -5./2.*logd[IDX(i,j,k)])/(dx[i]*dx[i])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(max(i-2,0),j,k)]-logd[IDX(min(i+2,in-1),j,k)])
          				  -2./3.*(logd[IDX(max(i-1,0),j,k)]-logd[IDX(min(i+1,in-1),j,k)]))/dx[i],2)/4.;
#endif /* QP_4TH_ORDER */

          //visx[IDX(i,j,k)] = u[IDX(max(i-1,0),j,k)]-2.0*u[IDX(i,j,k)]+u[IDX(min(i+1,in-1),j,k)];

          if (rank > 1) {
#ifdef QP_1ST_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + (pow(d[IDX(i,max(j-1,0),k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,min(j+1,jn-1),k)],0.5))/pow(dy[j],2);
#endif

#ifdef QP_2ND_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + (logd[IDX(i,max(j-1,0),k)]-2.0*logd[IDX(i,j,k)]+logd[IDX(i,min(j+1,jn-1),k)])/(dy[j]*dy[j])/2.;
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(i,min(j+1,jn-1),k)]-logd[IDX(i,max(j-1,0),k)])/(dy[j]*2.),2)/4.;
#endif

	    //4th order differentiation
#ifdef QP_4TH_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,max(j-2,0),k)]+logd[IDX(i,min(j+2,jn-1),k)])
          				  +4./3.*(logd[IDX(i,max(j-1,0),k)]+logd[IDX(i,min(j+1,jn-1),k)])
          				  -5./2.*logd[IDX(i,j,k)])/(dy[j]*dy[j])/2.;
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,max(j-2,0),k)]-logd[IDX(i,min(j+2,jn-1),k)])
          				  -2./3.*(logd[IDX(i,max(j-1,0),k)]-logd[IDX(i,min(j+1,jn-1),k)]))/dy[j],2)/4.;
#endif /* QP_4TH_ORDER */

          }// end rank > 1

          if (rank > 2) {
#ifdef QP_1ST_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + (pow(d[IDX(i,j,max(k-1,0))],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,j,min(k+1,kn-1))],0.5))/pow(dz[k],2);
#endif
		
#ifdef QP_2ND_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + (logd[IDX(i,j,max(k-1,0))]-2.0*logd[IDX(i,j,k)]+logd[IDX(i,j,min(k+1,kn-1))])/(dz[k]*dz[k])/2.;
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(i,j,min(k+1,kn-1))]-logd[IDX(i,j,max(k-1,0))])/(dz[k]*2.),2)/4.;
#endif
	    
          //4th order differentiation
#ifdef QP_4TH_ORDER
	    p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,j,max(k-2,0))]+logd[IDX(i,j,min(k+2,kn-1))])
          				  +4./3.*(logd[IDX(i,j,max(k-1,0))]+logd[IDX(i,j,min(k+1,kn-1))])
          				  -5./2.*logd[IDX(i,j,k)])/(dz[k]*dz[k])/2.;
	    p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,j,max(k-2,0))]-logd[IDX(i,j,min(k+2,kn-1))])
          				  -2./3.*(logd[IDX(i,j,max(k-1,0))]-logd[IDX(i,j,min(k+1,kn-1))]))/dz[k],2)/4.;
#endif /* QP_4TH_ORDER */

          }// end rank > 2

          p[IDX(i,j,k)] = p[IDX(i,j,k)]*lapcoef;
          
          e[IDX(i,j,k)] = 1e3;          

      }// end loop over i
    } // end loop over j
  } // end loop over k 

  /* 1) Substep 1 -- pressure and gravity terms */
      /* Update velocities */

  for (k = ksm1; k <= kep2; k++) {
    for (j = jsm1; j <= jep2; j++) {
      for (i = is-1; i <= ie+2; i++) {
      		
	deltav =  dt*(p[IDX(i,j,k)]-p[IDX(i-1,j,k)])/dx[i];
	u[IDX(i,j,k)] = u[IDX(i,j,k)] + deltav;

	if (rank > 1) {
	  deltav =  dt*(p[IDX(i,j,k)]-p[IDX(i,j-1,k)])/dy[j];
      	  v[IDX(i,j,k)] = v[IDX(i,j,k)] + deltav;
	}// end rank >1

	if (rank > 2) {
	  deltav =  dt*(p[IDX(i,j,k)]-p[IDX(i,j,k-1)])/dz[k];
      	  w[IDX(i,j,k)] = w[IDX(i,j,k)] + deltav;
        }// end rank>2

	/* Update velocities with acceleration */

        if(gravity ==1 ){

    	  u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*gr_xacc[IDX(i,j,k)];

	  if (rank > 1)
	    v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*gr_yacc[IDX(i,j,k)];

	  if (rank > 2)
	    w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*gr_zacc[IDX(i,j,k)];
	}
	e[IDX(i,j,k)] = 1e3;

      }// end loop over i
    } // end: loop over j
  } // end: loop over k

  delete [] logd;

  return SUCCESS;
}
