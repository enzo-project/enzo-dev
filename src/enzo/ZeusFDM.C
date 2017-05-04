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
  //float *visx = new float[size];
  //float *visy = new float[size];
  //float *visz = new float[size];
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
        u[i],v[i],w[i],d[i],e[i], dx[0], dt);
      ENZO_FAIL("ZeusSolver: Negative Density! \n");
    } else{
      logd[i] = log(d[i]);
    }
  }

  /* Compute the quantum pressure */
  

 /* for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){

          //p[IDX(i,j,k)] = (pow(d[IDX(max(i-1,0),j,k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(min(i+1,in-1),j,k)],0.5))/pow(dx[i],2);

          //2nd order differentiation
          //p[IDX(i,j,k)] = (logd[IDX(max(i-1,0),j,k)]-2.0*logd[IDX(i,j,k)]+logd[IDX(min(i+1,in-1),j,k)])/(dx[i]*dx[i])/2.;
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(min(i+1,in-1),j,k)]-logd[IDX(max(i-1,0),j,k)])/(dx[i]*2.),2)/4.;


          //4th order differentiation
          p[IDX(i,j,k)] = (-1./12.*(logd[IDX(max(i-2,0),j,k)]+logd[IDX(min(i+2,in-1),j,k)])
          				  +4./3.*(logd[IDX(max(i-1,0),j,k)]+logd[IDX(min(i+1,in-1),j,k)])
          				  -5./2.*logd[IDX(i,j,k)])/(dx[i]*dx[i])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(max(i-2,0),j,k)]-logd[IDX(min(i+2,in-1),j,k)])
          				  -2./3.*(logd[IDX(max(i-1,0),j,k)]-logd[IDX(min(i+1,in-1),j,k)]))/dx[i],2)/4.;

          //visx[IDX(i,j,k)] = u[IDX(max(i-1,0),j,k)]-2.0*u[IDX(i,j,k)]+u[IDX(min(i+1,in-1),j,k)];

          if (rank > 1) {
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + (pow(d[IDX(i,max(j-1,0),k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,min(j+1,jn-1),k)],0.5))/pow(dy[j],2);
		 
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + (logd[IDX(i,max(j-1,0),k)]-2.0*logd[IDX(i,j,k)]+logd[IDX(i,min(j+1,jn-1),k)])/(dy[j]*dy[j])/2.;
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(i,min(j+1,jn-1),k)]-logd[IDX(i,max(j-1,0),k)])/(dy[j]*2.),2)/4.;

          	//4th order differentiation
          p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,max(j-2,0),k)]+logd[IDX(i,min(j+2,jn-1),k)])
          				  +4./3.*(logd[IDX(i,max(j-1,0),k)]+logd[IDX(i,min(j+1,jn-1),k)])
          				  -5./2.*logd[IDX(i,j,k)])/(dy[j]*dy[j])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,max(j-2,0),k)]-logd[IDX(i,min(j+2,jn-1),k)])
          				  -2./3.*(logd[IDX(i,max(j-1,0),k)]-logd[IDX(i,min(j+1,jn-1),k)]))/dy[j],2)/4.;

          //visy[IDX(i,j,k)] = v[IDX(max(i-1,0),j,k)]-2.0*v[IDX(i,j,k)]+v[IDX(min(i+1,in-1),j,k)];

          //visx[IDX(i,j,k)] += u[IDX(i,max(j-1,0),k)]-2.0*u[IDX(i,j,k)]+u[IDX(i,min(j+1,jn-1),k)];
          //visy[IDX(i,j,k)] += v[IDX(i,max(j-1,0),k)]-2.0*v[IDX(i,j,k)]+v[IDX(i,min(j+1,jn-1),k)];



          }// end rank > 1
          if (rank > 2) {
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + (pow(d[IDX(i,j,max(k-1,0))],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,j,min(k+1,kn-1))],0.5))/pow(dz[k],2);
		
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + (logd[IDX(i,j,max(k-1,0))]-2.0*logd[IDX(i,j,k)]+logd[IDX(i,j,min(k+1,kn-1))])/(dz[k]*dz[k])/2.;
          //p[IDX(i,j,k)] = p[IDX(i,j,k)] + pow((logd[IDX(i,j,min(k+1,kn-1))]-logd[IDX(i,j,max(k-1,0))])/(dz[k]*2.),2)/4.;

          //4th order differentiation
          p[IDX(i,j,k)] = p[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,j,max(k-2,0))]+logd[IDX(i,j,min(k+2,kn-1))])
          				  +4./3.*(logd[IDX(i,j,max(k-1,0))]+logd[IDX(i,j,min(k+1,kn-1))])
          				  -5./2.*logd[IDX(i,j,k)])/(dz[k]*dz[k])/2.;
          p[IDX(i,j,k)] = p[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,j,max(k-2,0))]-logd[IDX(i,j,min(k+2,kn-1))])
          				  -2./3.*(logd[IDX(i,j,max(k-1,0))]-logd[IDX(i,j,min(k+1,kn-1))]))/dz[k],2)/4.;

          //visz[IDX(i,j,k)] = w[IDX(max(i-1,0),j,k)]-2.0*w[IDX(i,j,k)]+w[IDX(min(i+1,in-1),j,k)];

          //visz[IDX(i,j,k)] += w[IDX(i,max(j-1,0),k)]-2.0*w[IDX(i,j,k)]+w[IDX(i,min(j+1,jn-1),k)];

          //visx[IDX(i,j,k)] += u[IDX(i,j,max(k-1,0))]-2.0*u[IDX(i,j,k)]+u[IDX(i,j,min(k+1,kn-1))];
          //visy[IDX(i,j,k)] += v[IDX(i,j,max(k-1,0))]-2.0*v[IDX(i,j,k)]+v[IDX(i,j,min(k+1,kn-1))];
          //visz[IDX(i,j,k)] += w[IDX(i,j,max(k-1,0))]-2.0*w[IDX(i,j,k)]+w[IDX(i,j,min(k+1,kn-1))];


          }// end rank > 2

          //p[IDX(i,j,k)] = -1.0*p[IDX(i,j,k)]/pow(d[IDX(i,j,k)],0.5)*lapcoef;
          
          p[IDX(i,j,k)] = -1.0*p[IDX(i,j,k)]*lapcoef;
          
          e[IDX(i,j,k)] = 1e3;          
                  }// end loop over i
    } // end loop over j
  } // end loop over k

  /* 1) Substep 1 -- pressure and gravity terms */
      /* Update velocities */

  for (k = ksm1; k <= kep2; k++) {
    for (j = jsm1; j <= jep2; j++) {
      for (i = is-1; i <= ie+2; i++) {
      		
          /*deltav =  dt*(p[IDX(i-1,j,k)]-p[IDX(i,j,k)])/dx[i];// + visx[IDX(i,j,k)]* 1e-1;
		      u[IDX(i,j,k)] = u[IDX(i,j,k)] + deltav;

      if (rank > 1) {
      		deltav =  dt*(p[IDX(i,j-1,k)]-p[IDX(i,j,k)])/dy[j];// + visy[IDX(i,j,k)]* 1e-1;
      	  v[IDX(i,j,k)] = v[IDX(i,j,k)] + deltav;
          }// end rank >1

      if (rank > 2) {
      		deltav =  dt*(p[IDX(i,j,k-1)]-p[IDX(i,j,k)])/dz[k];// + visz[IDX(i,j,k)]* 1e-1;
      	  w[IDX(i,j,k)] = w[IDX(i,j,k)] + deltav;
        }// end rank>2

      /* Update velocities with acceleration */

        //if(gravity ==1 ){

    	  u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*gr_xacc[IDX(i,j,k)];

	     if (rank > 1)
  	       v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*gr_yacc[IDX(i,j,k)];

	     if (rank > 2)
	         w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*gr_zacc[IDX(i,j,k)];
       //}

      }// end loop over i
    } // end: loop over j
  } // end: loop over k
 /* 2) Substep 2 -- artificial viscosity */

  //      nstep = 1;
  nstep = 5;
  dt1 = dt/float(nstep);

  for (n = 0; n < nstep; n++) {
    for (k = ksm2; k < kn; k++) {
      for (j = jsm2; j < jn; j++) {

	/* a) Quadratic viscosity */

	for (i = is-2; i <= ie+2; i++) {
	  if ((u[IDX(i+1,j,k)]-u[IDX(i,j,k)]) < 0.0){
	    q[i] = C2*d[IDX(i,j,k)]*pow(u[IDX(i+1,j,k)]-u[IDX(i,j,k)],2.0);}
	  else{
	    q[i] = 0.0;}
	}

	/* b) linear viscosity */

	if (C1 != 0.0) {
	    for (i = is-2; i <= ie+2; i++){
	      q[i] = q[i] + C1*d[IDX(i,j,k)]*(u[IDX(i+1,j,k)]-u[IDX(i,j,k)])*
		    sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);}
	}

	q[is-3] = q[is-2];
	q[ie+3] = q[ie+2];

	/* update velocity1 and energy */

	for (i = is-2; i <= ie+3; i++){
	  u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt1*(q[i-1]-q[i])/
		(dx[i]*0.5*(d[IDX(i,j,k)]+d[IDX(i-1,j,k)]));}

      } // end loop over j
    } // end loop over k

    /* update velocity2 and energy */

    if (rank > 1) {
      for (k = ksm2; k < kn; k++) {
	for (i = is-2; i < in; i++) {

	  for (j = js-2; j <= je+2; j++) {
	    if ((v[IDX(i,j+1,k)]-v[IDX(i,j,k)]) < 0.0){
	      q[j] = C2*d[IDX(i,j,k)]*pow(v[IDX(i,j+1,k)]-v[IDX(i,j,k)], 2.0);
	    }
	    else{
	      q[j] = 0.0;
	    }
	  }

	  if (C1 != 0.0) {
	      for (j = js-2; j <= je+2; j++){
		q[j] = q[j] + C1*d[IDX(i,j,k)]*(v[IDX(i,j+1,k)]-v[IDX(i,j,k)])*
		      sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    }
	  }

	  q[js-3] = q[js-2];
	  q[je+3] = q[je+2];

	  for (j = js-2; j <= je+3; j++){
	    v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt1*(q[j-1]-q[j])/
 		            (dy[j]*0.5*(d[IDX(i,j,k)]+d[IDX(i,j-1,k)]));}

	} // end: loop over i
      } // end: loop over k
    } // end: if (rank > 1)

    /*  update velocity3 and energy */

    if (rank > 2) {
      for (j = jsm2; j < jn; j++) {
	for (i = is-2; i < in; i++) {

	  for (k = ks-2; k <= ke+2; k++) {
	    if ((w[IDX(i,j,k+1)]-w[IDX(i,j,k)]) < 0.0){
	      q[k] = C2*d[IDX(i,j,k)]*pow(w[IDX(i,j,k+1)]-w[IDX(i,j,k)], 2.0);
	    }
	    else{
	      q[k] = 0.0;
	    }
	  }

	  if (C1 != 0.0) {
	      for (k = ks-2; k <= ke+2; k++){
		q[k] = q[k] + C1*d[IDX(i,j,k)]*(w[IDX(i,j,k+1)]-w[IDX(i,j,k)])*
		                  sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    }
	  }

	  q[ks-3] = q[ks-2];
	  q[ke+3] = q[ke+2];


	  for (k = ks-2; k <= ke+3; k++){
	    w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt1*(q[k-1]-q[k])/
		            (dz[k]*0.5*(d[IDX(i,j,k)]+d[IDX(i,j,k-1)]));}

	} // end loop over i
      } // end loop over j
    } // end: if (rank > 2)

  } // end: loop over nstep (end of artificial viscosity)

  /* 3) Add MinimumPressureSupport */

  if (minsupecoef>0){
    // compute pressure
  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      for (i = 0; i < in; i++){
    e1 = minsupecoef*d[IDX(i,j,k)];
    gamma1 = 5./3.;
    p[IDX(i,j,k)] = (gamma1-1.0)*d[IDX(i,j,k)]*e1;
                              }// end loop over i
                            }// end loop over j
                          }// end loop over k

  for (k = ksm2; k <= kep2; k++) {
    for (j = jsm2; j <= jep2; j++) {
      for (i = is-1; i <= ie+2; i++) {
          
          deltav =  dt*(p[IDX(i-1,j,k)]-p[IDX(i,j,k)])/(dx[i]*0.5*(d[IDX(i-1,j,k)]+d[IDX(i,j,k)]));
          u[IDX(i,j,k)] = u[IDX(i,j,k)] + deltav;

      if (rank > 1) {
          deltav = dt*(p[IDX(i,j-1,k)]-p[IDX(i,j,k)])/(dy[j]*0.5*(d[IDX(i,j-1,k)]+d[IDX(i,j,k)]));
          v[IDX(i,j,k)] = v[IDX(i,j,k)] + deltav;
          }
      if (rank > 2) {
          deltav = dt*(p[IDX(i,j,k-1)]-p[IDX(i,j,k)])/(dz[k]*0.5*(d[IDX(i,j,k-1)]+d[IDX(i,j,k)]));
          w[IDX(i,j,k)] = w[IDX(i,j,k)] + deltav;
          }
        }// end loop over i
      }// end loop over j
    }// end loop over k

          }// end min pressure support




delete [] logd;

//delete [] visx;
//delete [] visy;
//delete [] visz;



  return SUCCESS;
}
