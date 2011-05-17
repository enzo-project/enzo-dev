/***********************************************************************
c
c  TRANSPORT TERM IN X-DIRECTION FOR ZEUS HYDRO (CARTESIAN ONLY)
c
c  written by: Greg Bryan (implemented from Stone & Norman, ApJS 80, 753)
c  date:       February, 1997
c  modified1: 
c
c  PURPOSE:
c
c  EXTERNALS:
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - gas specific energy field
c     i,j,kn  - dimensions of field arrays
c     rank    - dimension of problem (not currently used)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "Fluxes.h"
#include "fortran.def"

#define IDX(a,b,c) (( (c)*jn + (b) )*in + (a))

int vanlr1_zc(float *q, int idim, int jdim, int kdim, int is, int ie, 
	      int j, int k, float dx[], float dt, float u[], float qstar[]);
int vanlr1_fc(float *q, int idim, int jdim, int kdim, int is, int ie, 
	      int j, int k, float dx[], float dt, float u[], float qstar[]);


int Zeus_xTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dx[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int ncolor, int colnum[])
		    
{
  /*  Parameters */

  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Locals */

  int i, j, k, jm1, km1, i1, i2, n, idim, jdim, dim, offset, ic,
      fistart, fjstart, fiend, fjend, lface, rface;
  float dnew, q[ijk], div[ijk], f2[ijk], f3[ijk], f4[ijk], f5[ijk],
        dstar[ijk], estar[ijk], ustar[ijk], vstar[ijk], wstar[ijk],
        uavgi[ijk], uavgj[ijk], uavgk[ijk], df, ueff[ijk];
  float *colstar[MAX_COLOR];
  int ip1, im2, im3, ip2;

//=======================================================================

  /* Set the index for the in-scan (dim) off-scan (i,jdim) dimensions */

  dim = 0, idim = 1, jdim = 2;

  /* Allocate space for colstar. */

  for (ic = 0; ic < ncolor; ic++)
    colstar[ic] = new float[in];

// 1) Transport step - x direction

//  Compute the mass flux over the entire field

  for (k = 0; k < kn; k++) {
    for (j = 0; j < jn; j++) {

//        Interpolate density
      
      for (i = 0; i < in; i++)
	ueff[i] = u[IDX(i,j,k)];

      vanlr1_zc(d, in, jn, kn, is-1, ie+2, j, k, dx, dt, ueff, dstar);

//        Compute mass flux
      
      for (i = is-2; i <= ie+3; i++) {
	f1[IDX(i,j,k)] = dstar[i]*ueff[i];
	if (fabs(u[IDX(i,j,k)]) > 0.5*dx[i]/dt) {
	  printf("xt warning 1: u=%"GSYM",%"GSYM",%"GSYM"  dstar=%"GSYM",%"GSYM",%"GSYM"  i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
                 u[IDX(i-1,j,k)],u[IDX(i,j,k)],u[IDX(i+1,j,k)], 
		 dstar[i-1], dstar[i],dstar[i+1],i,j,k);
	  printf("  d=%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM",%"GSYM"   dt=%"GSYM"\n",
		 d[IDX(i-1,j,k)],d[IDX(i,j,k)],d[IDX(i+1,j,k)],
		 dx[i],dx[i+1],dt);
	}
      }

    } // end loop over j
  } // end loop over k

//  Update the quantities (loop backwards on j,k so j-1,k-1 are old values)

  for (k = kn-1; k >= 0; k--) {
    for (j = jn-1; j >= 0; j--) {
      jm1 = max(j-1, 0);
      km1 = max(k-1, 0);

      for (i = 0; i < in; i++)
	ueff[i] = u[IDX(i,j,k)];

//    Interpolate energy

      vanlr1_zc(e, in, jn, kn, is-1, ie+2, j, k, dx, dt, ueff, estar);

//    Compute color flux (assuming it is density-like) and advect
//           color variables (zone-centered)

      for (ic = 0; ic < ncolor; ic++) {

	vanlr1_zc(BaryonField[colnum[ic]], in, jn, kn, is-1, ie+2, j, k, dx, dt, ueff, colstar[ic]);

	for (i = is-2; i <= ie+3; i++)
	  colstar[ic][i] *= ueff[i];
	    
	for (i = is-2; i <= ie+2; i++) 
	  BaryonField[colnum[ic]][IDX(i,j,k)] += dt*(colstar[ic][i] - colstar[ic][i+1])/dx[i];
      }

//    Compute energy flux

      for (i = is-2; i <= ie+3; i++)
	f5[i] = estar[i]*f1[IDX(i,j,k)];

//    Make appropriately-averaged quanitities for advection

      for (i=is-3; i <= ie+2; i++) 
	uavgi[i] = 0.5*(ueff[i] + ueff[i+1]); // yes, this is right
      for (i=is-3; i <= ie+3; i++) {
	uavgj[i] = 0.5*(u[IDX(i,j,k)] + u[IDX(i,jm1,k)]);
	uavgk[i] = 0.5*(u[IDX(i,j,k)] + u[IDX(i,j,km1)]);
      }

//    Interpolate velocities

      vanlr1_fc(u, in, jn, kn, is-1, ie+1, j, k, dx, dt, uavgi, ustar);
      if (rank > 1) vanlr1_zc(v, in, jn, kn, is-1,ie+2,j,k, dx, dt, uavgj, vstar);
      if (rank > 2) vanlr1_zc(w, in, jn, kn, is-1,ie+2,j,k, dx, dt, uavgk, wstar);

//    Compute momentum fluxes

      for (i = is-2; i <= ie+2; i++)
	f2[i] = ustar[i]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i+1,j,k)]);

      if (rank > 1)
	for (i = is-1; i <= ie+2; i++) {
	  f3[i] = vstar[i]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,jm1,k)]);
	  f4[i] = wstar[i]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,j,km1)]);
	}

//    Convert velocities to momenta

      for (i = is-1; i <= ie+2; i++) {
	u[IDX(i,j,k)] = u[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i-1,j,k)]);
	u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*(f2[i-1]-f2[i])/dx[i];
      }
      if (rank > 1) {
	for (i = is-1; i <= ie+1; i++){
	  v[IDX(i,j,k)] = v[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,jm1,k)]);
	  w[IDX(i,j,k)] = w[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,j,km1)]);

	  /* Update momentum fluxes */

	  v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*(f3[i]-f3[i+1])/dx[i];
	  w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*(f4[i]-f4[i+1])/dx[i];
	}
      }

//    Update mass and energy fluxes

      for (i = is-2; i <= ie+2; i++) {

	dnew = d[IDX(i,j,k)] + dt*(f1[IDX(i,j,k)] - f1[IDX(i+1,j,k)])/dx[i];

	if (dnew/d[IDX(i,j,k)] < 0.15)
	  printf("zeus_x d warning: d,dnew=%"GSYM",%"GSYM"  i,j,k=%"ISYM",%"ISYM",%"ISYM"  f1=%"GSYM",%"GSYM"\n",
		 d[IDX(i,j,k)],dnew,i,j,k,f1[IDX(i,j,k)],f1[IDX(i+1,j,k)]);

	e[IDX(i,j,k)] = (e[IDX(i,j,k)]*d[IDX(i,j,k)] + dt*(f5[i] - f5[i+1])/dx[i])/dnew;

	if (e[IDX(i,j,k)] <= 0.0 || dnew <= 0.0) {
	  ENZO_VFAIL("zeus_x negative e or d error: d,e,dnew,dt=%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n",d[IDX(i,j,k)],e[IDX(i,j,k)],dnew,dt)
	}

	d[IDX(i,j,k)] = dnew;
      }

//    Convert momenta back to velocities

      for (i = is-1; i <= ie+2; i++) {

	u[IDX(i,j,k)] = u[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i-1,j,k)]));

	if (fabs(u[IDX(i,j,k)]) > dx[i]/dt) {
	  printf("zeus_x ux error: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke=%"ISYM",%"ISYM",%"ISYM"  dx,dt=%"GSYM",%"GSYM"\n", i,j,k,ie,je,ke,dx[i],dt);
          printf("               : u,d,d=%"GSYM",%"GSYM"%"GSYM"\n", u[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i-1,j,k)]);
	  for (i1=1; i1 <= in; i1++)
	    printf("%"ISYM" d,u,v,e,w,d-1,f2,u*,f1,uav=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n", 
		   i1, d[IDX(i1,j,k)],u[IDX(i1,j,k)],v[IDX(i1,j,k)], e[IDX(i1,j,k)],w[IDX(i1,j,k)],d[IDX(i1-1,j,k)], f2[i1],ustar[i1],f1[IDX(i1,j,k)],uavgi[i1]);
	  ENZO_FAIL("Velocity too fast!\n");
	}
      }

//  Check this slice against the list of subgrids 
//  (all subgrid quantities are zero based)

      for (n=0; n < nsubgrids; n++) {

	fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] -  
	  GridGlobalStart[idim];
	fiend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
	  GridGlobalStart[idim];
	fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
	  GridGlobalStart[jdim];
	fjend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
	  GridGlobalStart[jdim];
        lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] -
          GridGlobalStart[dim];
        rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
          GridGlobalStart[dim] + 1;

	if (k >= fjstart && k <= fjend && j >= fistart && j <= fiend) {

	  offset = j-fistart + (k-fjstart)*(fiend - fistart + 1);
	  i1 = lface;
	  i2 = rface;
	  SubgridFluxes[n]->LeftFluxes[DensNum][0][offset]  = f1[IDX(i1,j,k)]*dt/dx[i1];
	  SubgridFluxes[n]->RightFluxes[DensNum][0][offset] = f1[IDX(i2,j,k)]*dt/dx[i2];
//	  SubgridFluxes[n]->LeftFluxes[TENum][0][offset]    = f5[i1]*dt;
//        SubgridFluxes[n]->RightFluxes[TENum][0][offset]   = f5[i2]*dt;
//	  SubgridFluxes[n]->LeftFluxes[Vel1Num][0][offset]  = f2[i1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel1Num][0][offset] = f2[i2]*dt
//        if (rank > 1) {
//	    SubgridFluxes[n]->LeftFluxes[Vel2Num][0][offset]  = f3[i1]*dt;
//          SubgridFluxes[n]->RightFluxes[Vel2Num][0][offset] = f3[i2]*dt;
//        }
//        if (rank > 2) {
//	    SubgridFluxes[n]->LeftFluxes[Vel3Num][0][offset]  = f4[i1]*dt;
//          SubgridFluxes[n]->RightFluxes[Vel3Num][0][offset] = f4[i2]*dt;
//        }
	  for (ic=0; ic < ncolor; ic++) {
	    SubgridFluxes[n]->LeftFluxes[colnum[ic]][0][offset] = colstar[ic][i1]*dt;
	    SubgridFluxes[n]->RightFluxes[colnum[ic]][0][offset] = colstar[ic][i2]*dt;
	  }
	}
      } // end loop over subgrids

    } // next j line

  } // next k line

//  Convert v,w momenta back to velocities

  if (rank > 1) {
    for (k = 0; k < kn; k++) {
      for (j = 0; j < jn; j++) {
	jm1 = max(j-1, 0);
	km1 = max(k-1, 0);
	for (i = is-1; i <= ie+1; i++) {
	  v[IDX(i,j,k)] = v[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,jm1,k)]));
	  w[IDX(i,j,k)] = w[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,j,km1)]));
	  if (fabs(v[IDX(i,j,k)]) > dx[i]/dt) {
	    printf("zeus_x vx warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,jm1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		   i,j,k,ie,je,ke,jm1);
	    printf("zeus_x vx warning: v,d,d-1=%"GSYM",%"GSYM",%"GSYM"  dx,dt=%"GSYM",%"GSYM"\n", v[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,jm1,k)], dx[i], dt);
	  }
	  if (fabs(w[IDX(i,j,k)]) > dx[i]/dt) {
	    printf("zeus_x wx warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,km1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		   i,j,k,ie,je,ke,km1);
	    printf("zeus_x wx warning: v,d,d-1=%"GSYM",%"GSYM",%"GSYM"\n", w[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,j,km1)]);
	  }
	} // end: loop over i
      } // end: loop over j
    } // end: loop over k
  } // end: if (rank > 1)


  /* Cleanup */

  for (ic = 0; ic < ncolor; ic++)
    delete [] colstar[ic];

  return SUCCESS;

}
