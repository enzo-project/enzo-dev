/***********************************************************************
c
c  TRANSPORT TERM IN Z-DIRECTION FOR ZEUS HYDRO (CARTESIAN ONLY)
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

int vanlr3_zc(float *q, int idim, int jdim, int kdim, int ks, int ke, 
	      int i, int j, float dz[], float dt, float w[], float qstar[]);
int vanlr3_fc(float *q, int idim, int jdim, int kdim, int ks, int ke, 
	      int i, int j, float dz[], float dt, float w[], float qstar[]);


int Zeus_zTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dz[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int ncolor, int colnum[])
{
  /*  Parameters */

  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Locals */

  int i, j, k, im1, jm1, k1, k2, n, dim, idim, jdim, offset, ic,
      fistart, fjstart, fiend, fjend, lface, rface;
  float dnew, q[ijk], div[ijk], f2[ijk], f3[ijk], f4[ijk], f5[ijk],
        dstar[ijk], estar[ijk], ustar[ijk], vstar[ijk], wstar[ijk],
        uavgi[ijk], uavgj[ijk], uavgk[ijk], df, ueff[ijk];
  float *colstar[MAX_COLOR];
  int km1, kp1, km2, km3, kp2;

//=======================================================================

  /* Set the index for the in-scan (dim) off-scan (i,jdim) dimensions */

  dim = 2, idim = 0, jdim = 1;

  /* Allocate space for colstar. */

  for (ic = 0; ic < ncolor; ic++)
    colstar[ic] = new float[kn];

// 1) Transport step - x direction

//  Compute the mass flux over the entire field

  for (j = 0; j < jn; j++) {
    for (i = 0; i < in; i++) {

//    Interpolate density
      
      for (k = 0; k < kn; k++)
	ueff[k] = w[IDX(i,j,k)];

      vanlr3_zc(d, in, jn, kn, ks-1, ke+2, i, j, dz, dt, ueff, dstar);

//    Compute mass flux
      
      for (k = ks-2; k <= ke+3; k++) {
	km1 = max(k-1, ks-2);
	f1[IDX(i,j,k)] = dstar[k]*ueff[k];
	if (fabs(w[IDX(i,j,k)]) > 0.5*dz[k]/dt) {
	  printf("zt problem 1: w=%"GSYM",%"GSYM",%"GSYM"  dstar=%"GSYM",%"GSYM",%"GSYM"  i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
                 w[IDX(i,j-1,k)],w[IDX(i,j,k)],w[IDX(i,j+1,k)], 
		 dstar[k-1], dstar[k],dstar[k+1],i,j,k);
	  printf("  d=%"GSYM",%"GSYM",%"GSYM"  dz=%"GSYM",%"GSYM"   dt=%"GSYM"\n",
		 d[IDX(i,j,k-1)],d[IDX(i,j,k)],d[IDX(i,j,k+1)],
		 dz[k],dz[k+1],dt);
	}
      }

    } // end loop over i
  } // end loop over j

//  Update the quantities (loop backwards on j,k so j-1,k-1 are old values)

  for (j = jn-1; j >= 0; j--) {
    for (i = in-1; i >= 0; i--) {
      im1 = max(i-1, 0);
      jm1 = max(j-1, 0);

      for (k = 0; k < kn; k++)
	ueff[k] = w[IDX(i,j,k)];

//    Interpolate energy

      vanlr3_zc(e, in, jn, kn, ks-1, ke+2, i, j, dz, dt, ueff, estar);

//    Compute color flux (assuming it is density-like) and advect
//           color variables (zone-centered)

      for (ic = 0; ic < ncolor; ic++) {

	vanlr3_zc(BaryonField[colnum[ic]], in, jn, kn, ks-1, ke+2, i, j, dz, dt, ueff, colstar[ic]);

	for (k = ks-2; k <= ke+3; k++)
	  colstar[ic][k] *= ueff[k];
	    
	for (k = ks-2; k <= ke+2; k++) 
	  BaryonField[colnum[ic]][IDX(i,j,k)] += dt*(colstar[ic][k] - colstar[ic][k+1])/dz[k];
      }

//    Compute energy flux

      for (k = ks-2; k <= ke+3; k++)
	f5[k] = estar[k]*f1[IDX(i,j,k)];

//    Make appropriately-averaged quanitities for advection

      for (k=ks-3; k <= ke+3; k++) {
	uavgi[k] = 0.5*(w[IDX(i,j,k)] + w[IDX(im1,j,k)]);
	uavgj[k] = 0.5*(w[IDX(i,j,k)] + w[IDX(i,jm1,k)]);
      }
      for (k=ks-3; k <= ke+2; k++) 
	uavgk[k] = 0.5*(ueff[k] + ueff[k+1]); // yes, this is right

//    Interpolate velocities

      vanlr3_zc(u, in, jn, kn, ks-1, ke+2, i, j, dz, dt, uavgi, ustar);
      vanlr3_zc(v, in, jn, kn, ks-1, ke+2, i, j, dz, dt, uavgj, vstar);
      vanlr3_fc(w, in, jn, kn, ks-1, ke+1, i, j, dz, dt, uavgk, wstar);

//    Compute momentum fluxes

      for (k = ks-1; k <= ke+2; k++) {
	f2[k] = ustar[k]*0.5*(f1[IDX(i,j,k)] + f1[IDX(im1,j,k)]);
	f3[k] = vstar[k]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,jm1,k)]);
      }

      for (k = ks-2; k <= ke+2; k++)
	f4[k] = wstar[k]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,j,k+1)]);

//    Convert velocities to momenta

      for (k = ks-1; k <= ke+2; k++) {
	w[IDX(i,j,k)] = w[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,j,k-1)]);
	w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*(f4[k-1]-f4[k])/dz[k];
      }
      for (k = ks-1; k <= ke+1; k++){
	u[IDX(i,j,k)] = u[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(im1,j,k)]);
	v[IDX(i,j,k)] = v[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,jm1,k)]);

//    Update momentum fluxes

	u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*(f2[k]-f2[k+1])/dz[k];
	v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*(f3[k]-f3[k+1])/dz[k];
      }

//    Update mass and energy fluxes

      for (k = ks-2; k <= ke+2; k++) {

	dnew = d[IDX(i,j,k)] + dt*(f1[IDX(i,j,k)] - f1[IDX(i,j,k+1)])/dz[k];

	e[IDX(i,j,k)] = (e[IDX(i,j,k)]*d[IDX(i,j,k)] + dt*(f5[k] - f5[k+1])/dz[k])/dnew;

	if (e[IDX(i,j,k)] <= 0.0 || dnew <= 0.0) {
	  ENZO_VFAIL("zeus_z negative e or d error: d,e,dnew,dt=%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n",d[IDX(i,j,k)],e[IDX(i,j,k)],dnew,dt)
	}

	d[IDX(i,j,k)] = dnew;
      }

//    Convert momenta back to velocities

      for (k = ks-1; k <= ke+2; k++) {

	w[IDX(i,j,k)] = w[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,j,k-1)]));

	if (fabs(w[IDX(i,j,k)]) > dz[k]/dt) {
	  printf("zeus_z uy error: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke=%"ISYM",%"ISYM",%"ISYM"  dz,dt=%"GSYM",%"GSYM"\n", i,j,k,ie,je,ke,dz[k],dt);
          printf("               : w,d,d=%"GSYM",%"GSYM"%"GSYM"\n", w[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,j,k-1)]);
	  for (k1=0; k1 < kn; k1++)
	    printf("%"ISYM" d,u,v,e,w,d-1,f2,u*,f1,uav=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n", 
		   k1, d[IDX(i,j,k1)],u[IDX(i,j,k1)],v[IDX(i,j,k1)], e[IDX(i,j,k1)],w[IDX(i,j,k1)],d[IDX(i,j,k1-1)], f2[k1],ustar[k1],f1[IDX(i,j,k1)],uavgi[k1]);
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

	if (j >= fjstart && j <= fjend && i >= fistart && i <= fiend) {
	  offset = i-fistart + (j-fjstart)*(fiend - fistart + 1);
	  k1 = lface;
	  k2 = rface;
	  SubgridFluxes[n]->LeftFluxes[DensNum][2][offset]  = f1[IDX(i,j,k1)]*dt/dz[k1];
	  SubgridFluxes[n]->RightFluxes[DensNum][2][offset] = f1[IDX(i,j,k2)]*dt/dz[k2];
//	  SubgridFluxes[n]->LeftFluxes[TENum][2][offset]    = f5[k1]*dt;
//        SubgridFluxes[n]->RightFluxes[TENum][2][offset]   = f5[k2]*dt;
//	  SubgridFluxes[n]->LeftFluxes[Vel1Num][2][offset]  = f2[k1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel1Num][2][offset] = f2[k2]*dt
//	  SubgridFluxes[n]->LeftFluxes[Vel2Num][2][offset]  = f3[k1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel2Num][2][offset] = f3[k2]*dt;
//	  SubgridFluxes[n]->LeftFluxes[Vel3Num][2][offset]  = f4[k1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel3Num][2][offset] = f4[k2]*dt;
	  for (ic=0; ic < ncolor; ic++) {
	    SubgridFluxes[n]->LeftFluxes[colnum[ic]][2][offset] = colstar[ic][k1]*dt;
	    SubgridFluxes[n]->RightFluxes[colnum[ic]][2][offset] = colstar[ic][k2]*dt;
	  }
	}
      } // end loop over subgrids

    } // next i line

  } // next j line

//  Convert v,w momenta back to velocities

  for (j = 0; j < jn; j++) {
    for (i = 0; i < in; i++) {
      im1 = max(i-1, 0);
      jm1 = max(j-1, 0);
      for (k = ks-1; k <= ke+1; k++) {
	u[IDX(i,j,k)] = u[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(im1,j,k)]));
	v[IDX(i,j,k)] = v[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,jm1,k)]));
	if (fabs(u[IDX(i,j,k)]) > dz[k]/dt) {
	  printf("zeus_z ux warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,im1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		 i,j,k,ie,je,ke,im1);
	  printf("zeus_z ux warning: u,d,d-1=%"GSYM",%"GSYM",%"GSYM"\n", u[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(im1,j,k)]);
	}
	if (fabs(v[IDX(i,j,k)]) > dz[k]/dt) {

	  printf("zeus_z wx warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,jm1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		 i,j,k,ie,je,ke,jm1);
	  printf("zeus_z wx warning: v,d,d-1=%"GSYM",%"GSYM",%"GSYM"\n", v[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,jm1,k)]);
	}
      } // end: loop over k
    } // end: loop over i
  } // end: loop over j

  /* Cleanup */

  for (ic = 0; ic < ncolor; ic++)
    delete [] colstar[ic];

  return SUCCESS;

}
