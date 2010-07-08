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

int ZeusSource(float *d, float *e, float *u, float *v, float *w, float *p, 
	       int in, int jn, int kn, int rank, int igamfield,
	       int is, int ie, int js, int je, int ks, int ke, 
	       float C1, float C2, int ipresfree,
	       float *gamma, float dt, float pmin, float dx[], float dy[], float dz[],
	       int gravity, float *gr_xacc, float *gr_yacc, float *gr_zacc, 
	       int bottom, float minsupecoef)
{
  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Local declarations */

  int i, j, k, jsm1, ksm1, jep1, kep1, jsm2, ksm2, jep2, kep2, n, nstep;
  float  alpha, q[ijk], div[ijk], deltav, deltavmax, e1, gamma1, dt1;

  /* ======================================================================= */

  /* Compute varients on start indexes */

  jsm1 = max(js-1, 0);
  ksm1 = max(ks-1, 0);
  jep1 = min(je+1, jn-1);
  kep1 = min(ke+1, kn-1);
  jsm2 = max(js-2, 0);
  ksm2 = max(ks-2, 0);
  jep2 = min(je+2, jn-1);
  kep2 = min(ke+2, kn-1);

  gamma1 = gamma[0];  // if gamma is a scalar

  /* Compute the pressure */

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {
      if (ipresfree == 1) {

	/* Pressurefree - set pressure field to zero */

	for (i = 0; i < in; i++)
	  p[IDX(i,j,k)] = 0;

      } else if (igamfield == 1) {

	/* Compute pressure with variable gamma */

	for (i = 0; i < in; i++) {
	  e1 = e[IDX(i,j,k)];
	  e1 = max(e[IDX(i,j,k)], minsupecoef*d[IDX(i,j,k)]);
	  p[IDX(i,j,k)] = max((gamma[IDX(i,j,k)]-1.0)*d[IDX(i,j,k)]*e1, pmin);
	}
      } else {

	/* Compute pressure with constant gamma */

	for (i = 0; i < in; i++) {
	  e1 = e[IDX(i,j,k)];
	  e1 = max(e[IDX(i,j,k)], minsupecoef*d[IDX(i,j,k)]);
	  p[IDX(i,j,k)] = max((gamma1-1.0)*d[IDX(i,j,k)]*e1, pmin);

	  if (e[IDX(i,j,k)] <= 0.0 || d[IDX(i,j,k)] <= 0.0) {
	    printf("%"ISYM"\n", IDX(i,j,k));
	    fprintf(stderr, "zeus_source1: e,d=%"GSYM",%"GSYM" i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
		    e[IDX(i,j,k)],d[IDX(i,j,k)],i,j,k);
	    ENZO_FAIL("Negative energy or density!\n");
	  }
	}
      }


    } // end loop over j
  } // end loop over k

  /* 1) Substep 1 -- pressure and gravity terms */

  for (k = ksm2; k < kn; k++) {
    for (j = jsm2; j < jn; j++) {

      /* Update velocities with compression term
	 (limit increase to preventy pressure driven instability) */

      for (i = is-2; i <= ie+3; i++) {
	deltav =         dt*(p[IDX(i-1,j,k)]-p[IDX(i,j,k)])/ 
	         (dx[i]*0.5*(d[IDX(i-1,j,k)]+d[IDX(i,j,k)]));
	u[IDX(i,j,k)] = u[IDX(i,j,k)] + deltav;
      }

      if (rank > 1) {
	for (i = is-2; i <= ie+3; i++) {
	  deltav =         dt*(p[IDX(i,j-1,k)]-p[IDX(i,j,k)])/
	           (dy[j]*0.5*(d[IDX(i,j-1,k)]+d[IDX(i,j,k)]));
	  v[IDX(i,j,k)] = v[IDX(i,j,k)] + deltav;
	}
      }
      if (rank > 2) {
	for (i = is-2; i <= ie+3; i++) {
	  deltav =         dt*(p[IDX(i,j,k-1)]-p[IDX(i,j,k)])/
            	   (dz[k]*0.5*(d[IDX(i,j,k-1)]+d[IDX(i,j,k)]));
	  w[IDX(i,j,k)] = w[IDX(i,j,k)] + deltav;
	}
      }

      /* Update velocities with acceleration */

      if (gravity == 1) {

	for (i = is-2; i <= ie+3; i++)
	  u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*gr_xacc[IDX(i,j,k)];

	if (rank > 1)
	  for (i = is-2; i <= ie+3; i++)
	    v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*gr_yacc[IDX(i,j,k)];

	if (rank > 2)
	  for (i = is-2; i <= ie+3; i++)
	    w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*gr_zacc[IDX(i,j,k)];

      }
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
	  if ((u[IDX(i+1,j,k)]-u[IDX(i,j,k)]) < 0.0)
	    q[i] = C2*d[IDX(i,j,k)]*pow(u[IDX(i+1,j,k)]-u[IDX(i,j,k)],2.0);
	  else
	    q[i] = 0.0;
	}

	/* b) linear viscosity */

	if (C1 != 0.0) {
	  if (igamfield == 1) {
	    for (i = is-2; i <= ie+2; i++)
	      q[i] = q[i] + C1*d[IDX(i,j,k)]*(u[IDX(i+1,j,k)]-u[IDX(i,j,k)])*
		    sqrt(gamma[IDX(i,j,k)]*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	  } else {
	    for (i = is-2; i <= ie+2; i++)
	      q[i] = q[i] + C1*d[IDX(i,j,k)]*(u[IDX(i+1,j,k)]-u[IDX(i,j,k)])*
		    sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	  }
	}

	q[is-3] = q[is-2];
	q[ie+3] = q[ie+2];

	/* update velocity1 and energy */

	for (i = is-2; i <= ie+2; i++)
	  e[IDX(i,j,k)] = e[IDX(i,j,k)] + dt1*q[i]/d[IDX(i,j,k)]*
		(u[IDX(i,j,k)]-u[IDX(i+1,j,k)])/dx[i];

	for (i = is-2; i <= ie+3; i++)
	  u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt1*(q[i-1]-q[i])/
		(dx[i]*0.5*(d[IDX(i,j,k)]+d[IDX(i-1,j,k)]));

      } // end loop over j
    } // end loop over k

    /* update velocity2 and energy */

    if (rank > 1) {
      for (k = ksm2; k < kn; k++) {
	for (i = is-2; i < in; i++) {

	  for (j = js-2; j <= je+2; j++) {
	    if ((v[IDX(i,j+1,k)]-v[IDX(i,j,k)]) < 0.0)
	      q[j] = C2*d[IDX(i,j,k)]*pow(v[IDX(i,j+1,k)]-v[IDX(i,j,k)], 2.0);
	    else
	      q[j] = 0.0;
	  }

	  if (C1 != 0.0) {
	    if (igamfield == 1) {
	      for (j = js-2; j <= je+2; j++)
		q[j] = q[j] + C1*d[IDX(i,j,k)]*(v[IDX(i,j+1,k)]-v[IDX(i,j,k)])*
		      sqrt(gamma[IDX(i,j,k)]*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    } else {
	      for (j = js-2; j <= je+2; j++)
		q[j] = q[j] + C1*d[IDX(i,j,k)]*(v[IDX(i,j+1,k)]-v[IDX(i,j,k)])*
		      sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    }
	  }

	  q[js-3] = q[js-2];
	  q[je+3] = q[je+2];

	  for (j = js-2; j <= je+2; j++)
	    e[IDX(i,j,k)] = e[IDX(i,j,k)] + dt1*q[j]/d[IDX(i,j,k)]*
 		            (v[IDX(i,j,k)]-v[IDX(i,j+1,k)])/dy[j];

	  for (j = js-2; j <= je+3; j++)
	    v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt1*(q[j-1]-q[j])/
 		            (dy[j]*0.5*(d[IDX(i,j,k)]+d[IDX(i,j-1,k)]));

	} // end: loop over i
      } // end: loop over k
    } // end: if (rank > 1)

    /*  update velocity3 and energy */

    if (rank > 2) {
      for (j = jsm2; j < jn; j++) {
	for (i = is-2; i < in; i++) {

	  for (k = ks-2; k <= ke+2; k++) {
	    if ((w[IDX(i,j,k+1)]-w[IDX(i,j,k)]) < 0.0)
	      q[k] = C2*d[IDX(i,j,k)]*pow(w[IDX(i,j,k+1)]-w[IDX(i,j,k)], 2.0);
	    else
	      q[k] = 0.0;
	  }

	  if (C1 != 0.0) {
	    if (igamfield == 1) {
	      for (k = ks-2; k <= ke+2; k++) 
		q[k] = q[k] + C1*d[IDX(i,j,k)]*(w[IDX(i,j,k+1)]-w[IDX(i,j,k)])*
		                   sqrt(gamma[IDX(i,j,k)]*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    } else {
	      for (k = ks-2; k <= ke+2; k++)
		q[k] = q[k] + C1*d[IDX(i,j,k)]*(w[IDX(i,j,k+1)]-w[IDX(i,j,k)])*
		                  sqrt(gamma1*p[IDX(i,j,k)]/d[IDX(i,j,k)]);
	    }
	  }

	  q[ks-3] = q[ks-2];
	  q[ke+3] = q[ke+2];

	  for (k = ks-2; k <= ke+2; k++)
	    e[IDX(i,j,k)] = e[IDX(i,j,k)] + dt1*q[k]/d[IDX(i,j,k)]*
		            (w[IDX(i,j,k)]-w[IDX(i,j,k+1)])/dz[k];

	  for (k = ks-2; k <= ke+3; k++)
	    w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt1*(q[k-1]-q[k])/
		            (dz[k]*0.5*(d[IDX(i,j,k)]+d[IDX(i,j,k-1)]));

	} // end loop over i
      } // end loop over j
    } // end: if (rank > 2)

  } // end: loop over nstep (end of artificial viscosity)

  /*  3) Substep 3 -- compression term */

  for (k = ksm2; k <= kep1; k++) {
    for (j = jsm2; j <= jep1; j++) {

      /*  Compute the divergence (should use old u,v,w?) */

      for (i = is-2; i <= ie+1; i++)
	div[i] = (u[IDX(i+1,j,k)] - u[IDX(i,j,k)])/dx[i];

      if (rank > 1) {
	for (i = is-2; i <= ie+1; i++)
	  div[i] = div[i] + (v[IDX(i,j+1,k)] - v[IDX(i,j,k)])/dy[j];
      }

      if (rank > 2) {
	for (i = is-2; i <= ie+1; i++)
	  div[i] = div[i] + (w[IDX(i,j,k+1)] - w[IDX(i,j,k)])/dz[k];
      }

      /*  Update energy */

      for (i = is-2; i<= ie+1; i++) {
	if (igamfield == 1)
	  alpha = 0.5*dt*(gamma[IDX(i,j,k)] - 1.0)*div[i];
	else
	  alpha = 0.5*dt*(gamma1 - 1.0)*div[i];

	e[IDX(i,j,k)] = e[IDX(i,j,k)] * (1.0 - alpha)/(1.0 + alpha);

	if (e[IDX(i,j,k)] <= 0.0 || d[IDX(i,j,k)] <= 0.0) {
	  fprintf(stderr, "zeus_div error: e,alpha,div=%"GSYM",%"GSYM",%"GSYM"   i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
                  e[IDX(i,j,k)],alpha, div[i],i,j,k);
	  fprintf(stderr, "%"GSYM" %"GSYM" dt=%"GSYM" dx=%"GSYM"\n", d[IDX(i,j,k)],p[IDX(i,j,k)],dt,dx[i]);
	  fprintf(stderr, "p1=%"GSYM",%"GSYM",%"GSYM"\n", p[IDX(i+1,j,k)],p[IDX(i,j+1,k)],p[IDX(i,j,k+1)]);
	  fprintf(stderr, "p2=%"GSYM",%"GSYM",%"GSYM"\n", p[IDX(i-1,j,k)],p[IDX(i,j-1,k)],p[IDX(i,j,k-1)]);
	  fprintf(stderr, "d=%"GSYM",%"GSYM",%"GSYM"\n", d[IDX(i-1,j,k)],d[IDX(i,j-1,k)],d[IDX(i,j,k-1)]);
	  fprintf(stderr, "u=%"GSYM",%"GSYM",%"GSYM"\n", u[IDX(i-1,j,k)],u[IDX(i,j,k)],u[IDX(i+1,j,k)]);
	  fprintf(stderr, "v=%"GSYM",%"GSYM",%"GSYM"\n", v[IDX(i,j-1,k)],v[IDX(i,j,k)],v[IDX(i,j+1,k)]);
	  fprintf(stderr, "w=%"GSYM",%"GSYM",%"GSYM"\n", w[IDX(i,j,k-1)],w[IDX(i,j,k)],w[IDX(i,j,k+1)]);
	  fprintf(stderr, "grav_x=%"GSYM",%"GSYM",%"GSYM"\n", gr_xacc[IDX(i-1,j,k)],gr_xacc[IDX(i,j,k)], gr_xacc[IDX(i+1,j,k)]);
	  fprintf(stderr, "grav_y=%"GSYM",%"GSYM",%"GSYM"\n", gr_yacc[IDX(i,j-1,k)],gr_yacc[IDX(i,j,k)], gr_yacc[IDX(i,j+1,k)]);
	  fprintf(stderr, "grav_z=%"GSYM",%"GSYM",%"GSYM"\n", gr_zacc[IDX(i,j,k-1)],gr_zacc[IDX(i,j,k)], gr_zacc[IDX(i,j,k+1)]);
	  for (n=is; n <= ie; n++)
	    fprintf(stderr, "d,e,u,v,w=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n", d[IDX(n,j,k)],e[IDX(n,j,k)],u[IDX(n,j,k)],v[IDX(n,j,k)],w[IDX(n,j,k)]);
	  ENZO_FAIL("Negative energy or density!\n");
	}

      } // end: loop over i

    } // end: loop over j
  } // end: loop over k

  /*  Compute the pressure */

#ifdef UNUSED
  for (k = ksm2; k <= kep2; k++) {
    for (j = jsm2; j<= jep2; j++) {
      for (i = is-2; i <= ie+2; i++) {
	p[IDX(i,j,k)] = max((gamma-1.0)*d[IDX(i,j,k)]*e[IDX(i,j,k)], pmin);
	if (e[IDX(i,j,k)] <= 0.0 || d[IDX(i,j,k)] <= 0.0) {
	  ENZO_VFAIL("zeus_source2: e,d=%"GSYM",%"GSYM"  i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
		  e(i,j,j),d[IDX(i,j,k)],i,j,k)

	}
      }
    }
  }
#endif /* UNUSED */

  return SUCCESS;
}
