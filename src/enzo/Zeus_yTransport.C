/***********************************************************************
c
c  TRANSPORT TERM IN Y-DIRECTION FOR ZEUS HYDRO (CARTESIAN ONLY)
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

int vanlr2_zc(float *q, int idim, int jdim, int kdim, int js, int je, 
	      int i, int k, float dy[], float dt, float v[], float qstar[]);
int vanlr2_fc(float *q, int idim, int jdim, int kdim, int js, int je, 
	      int i, int k, float dy[], float dt, float v[], float qstar[]);


int Zeus_yTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dy[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int ncolor, int colnum[])
{
  /*  Parameters */

  int ijk = MAX_ANY_SINGLE_DIRECTION;

  /* Locals */

  int i, j, k, im1, km1, j1, j2, n, dim, idim, jdim, offset, ic,
      fistart, fjstart, fiend, fjend, lface, rface;
  float dnew, q[ijk], div[ijk], f2[ijk], f3[ijk], f4[ijk], f5[ijk],
        dstar[ijk], estar[ijk], ustar[ijk], vstar[ijk], wstar[ijk],
        uavgi[ijk], uavgj[ijk], uavgk[ijk], df, ueff[ijk];
  float *colstar[MAX_COLOR];
  int jp1, jm2, jm3, jp2;

//=======================================================================

  /* Set the index for the in-scan (dim) off-scan (i,jdim) dimensions */

  dim = 1, idim = 0, jdim = 2;

  /* Allocate space for colstar. */

  for (ic = 0; ic < ncolor; ic++)
    colstar[ic] = new float[jn];

// 1) Transport step - x direction

//  Compute the mass flux over the entire field

  for (k = 0; k < kn; k++) {
    for (i = 0; i < in; i++) {

//    Interpolate density
      
      for (j = 0; j < jn; j++)
	ueff[j] = v[IDX(i,j,k)];

      vanlr2_zc(d, in, jn, kn, js-1, je+2, i, k, dy, dt, ueff, dstar);

//    Compute mass flux
      
      for (j = js-2; j <= je+3; j++) {
	f1[IDX(i,j,k)] = dstar[j]*ueff[j];
	if (fabs(v[IDX(i,j,k)]) > 0.5*dy[j]/dt) {
	  printf("yt problem 1: v=%"GSYM",%"GSYM",%"GSYM"  dstar=%"GSYM",%"GSYM",%"GSYM"  i,j,k=%"ISYM",%"ISYM",%"ISYM"\n",
                 v[IDX(i,j-1,k)],v[IDX(i,j,k)],v[IDX(i,j+1,k)], 
		 dstar[j-1], dstar[j],dstar[j+1],i,j,k);
	  printf("  d=%"GSYM",%"GSYM",%"GSYM"  dy=%"GSYM",%"GSYM"   dt=%"GSYM"\n",
		 d[IDX(i,j-1,k)],d[IDX(i,j,k)],d[IDX(i,j+1,k)],
		 dy[j],dy[j+1],dt);
	}
      }

    } // end loop over i
  } // end loop over k

//  Update the quantities (loop backwards on j,k so j-1,k-1 are old values)

  for (k = kn-1; k >= 0; k--) {
    for (i = in-1; i >= 0; i--) {
      im1 = max(i-1, 0);
      km1 = max(k-1, 0);

      for (j = 0; j < jn; j++)
	ueff[j] = v[IDX(i,j,k)];

//    Interpolate energy

      vanlr2_zc(e, in, jn, kn, js-1, je+2, i, k, dy, dt, ueff, estar);

//    Compute color flux (assuming it is density-like) and advect
//           color variables (zone-centered)

      for (ic = 0; ic < ncolor; ic++) {

	vanlr2_zc(BaryonField[colnum[ic]], in, jn, kn, js-1, je+2, i, k, dy, dt, ueff, colstar[ic]);

	for (j = js-2; j <= je+3; j++)
	  colstar[ic][j] *= ueff[j];
	    
	for (j = js-2; j <= je+2; j++) 
	  BaryonField[colnum[ic]][IDX(i,j,k)] += dt*(colstar[ic][j] - colstar[ic][j+1])/dy[j];
      }

//    Compute energy flux

      for (j = js-2; j <= je+3; j++)
	f5[j] = estar[j]*f1[IDX(i,j,k)];

//    Make appropriately-averaged quanitities for advection

      for (j=js-3; j <= je+2; j++) 
	uavgj[j] = 0.5*(ueff[j] + ueff[j+1]); // yes, this is right
      for (j=js-3; j <= je+3; j++) {
	uavgi[j] = 0.5*(v[IDX(i,j,k)] + v[IDX(im1,j,k)]);
	uavgk[j] = 0.5*(v[IDX(i,j,k)] + v[IDX(i,j,km1)]);
      }

//    Interpolate velocities

      vanlr2_zc(u, in, jn, kn, js-1, je+2, i, k, dy, dt, uavgi, ustar);
      vanlr2_fc(v, in, jn, kn, js-1, je+1, i, k, dy, dt, uavgj, vstar);
      if (rank > 2) vanlr2_zc(w, in, jn, kn, js-1,je+2,i,k, dy, dt, uavgk, wstar);

//    Compute momentum fluxes

      for (j = js-2; j <= je+2; j++)
	f3[j] = vstar[j]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,j+1,k)]);

      for (j = js-1; j <= je+2; j++) {
	f2[j] = ustar[j]*0.5*(f1[IDX(i,j,k)] + f1[IDX(im1,j,k)]);
	f4[j] = wstar[j]*0.5*(f1[IDX(i,j,k)] + f1[IDX(i,j,km1)]);
      }

//    Convert velocities to momenta

      for (j = js-1; j <= je+2; j++) {
	v[IDX(i,j,k)] = v[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,j-1,k)]);
	v[IDX(i,j,k)] = v[IDX(i,j,k)] + dt*(f3[j-1]-f3[j])/dy[j];
      }
      for (j = js-1; j <= je+1; j++){
	u[IDX(i,j,k)] = u[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(im1,j,k)]);
	w[IDX(i,j,k)] = w[IDX(i,j,k)]*0.5*(d[IDX(i,j,k)] + d[IDX(i,j,km1)]);

//    Update momentum fluxes

	u[IDX(i,j,k)] = u[IDX(i,j,k)] + dt*(f2[j]-f2[j+1])/dy[j];
	w[IDX(i,j,k)] = w[IDX(i,j,k)] + dt*(f4[j]-f4[j+1])/dy[j];
      }

//    Update mass and energy fluxes

      for (j = js-2; j <= je+2; j++) {

	dnew = d[IDX(i,j,k)] + dt*(f1[IDX(i,j,k)] - f1[IDX(i,j+1,k)])/dy[j];

	e[IDX(i,j,k)] = (e[IDX(i,j,k)]*d[IDX(i,j,k)] + dt*(f5[j] - f5[j+1])/dy[j])/dnew;

	if (e[IDX(i,j,k)] <= 0.0 || dnew <= 0.0) {
	  ENZO_VFAIL("zeus_y negative e or d error: d,e,dnew,dt=%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n",d[IDX(i,j,k)],e[IDX(i,j,k)],dnew,dt)
	}

	d[IDX(i,j,k)] = dnew;
      }

//    Convert momenta back to velocities

      for (j = js-1; j <= je+2; j++) {

	v[IDX(i,j,k)] = v[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,j-1,k)]));

	if (fabs(v[IDX(i,j,k)]) > dy[j]/dt) {
	  printf("zeus_y uy error: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke=%"ISYM",%"ISYM",%"ISYM"  dy,dt=%"GSYM",%"GSYM"\n", i,j,k,ie,je,ke,dy[j],dt);
          printf("               : u,d,d=%"GSYM",%"GSYM"%"GSYM"\n", v[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,j-1,k)]);
	  for (j1=0; j1 < jn; j1++)
	    printf("%"ISYM" d,u,v,e,w,d-1,f2,u*,f1,uav=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"\n", 
		   j1, d[IDX(i,j1,k)],u[IDX(i,j1,k)],v[IDX(i,j1,k)], e[IDX(i,j1,k)],w[IDX(i,j1,k)],d[IDX(i,j1-1,k)], f2[j1],ustar[j1],f1[IDX(i,j1,k)],uavgi[j1]);
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

	if (i >= fistart && i <= fiend && k >= fjstart && k <= fjend) {

	  offset = i-fistart + (k-fjstart)*(fiend-fistart+1);
	  j1 = lface;
	  j2 = rface;
	  SubgridFluxes[n]->LeftFluxes[DensNum][1][offset]  = f1[IDX(i,j1,k)]*dt/dy[j1];
	  SubgridFluxes[n]->RightFluxes[DensNum][1][offset] = f1[IDX(i,j2,k)]*dt/dy[j2];
//	  SubgridFluxes[n]->LeftFluxes[TENum][1][offset]    = f5[j1]*dt;
//        SubgridFluxes[n]->RightFluxes[TENum][1][offset]   = f5[j2]*dt;
//	  SubgridFluxes[n]->LeftFluxes[Vel1Num][1][offset]  = f2[j1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel1Num][1][offset] = f2[j2]*dt
//	  SubgridFluxes[n]->LeftFluxes[Vel2Num][1][offset]  = f3[j1]*dt;
//        SubgridFluxes[n]->RightFluxes[Vel2Num][1][offset] = f3[j2]*dt;
//        if (rank > 2) {
//	    SubgridFluxes[n]->LeftFluxes[Vel3Num][1][offset]  = f4[j1]*dt;
//          SubgridFluxes[n]->RightFluxes[Vel3Num][1][offset] = f4[j2]*dt;
//        }
	  for (ic=0; ic < ncolor; ic++) {
	    SubgridFluxes[n]->LeftFluxes[colnum[ic]][1][offset] = colstar[ic][j1]*dt;
	    SubgridFluxes[n]->RightFluxes[colnum[ic]][1][offset] = colstar[ic][j2]*dt;
	  }
	}
      } // end loop over subgrids

    } // next i line

  } // next k line

//  Convert v,w momenta back to velocities

  for (k = 0; k < kn; k++) {
    for (i = 0; i < in; i++) {
      im1 = max(i-1, 0);
      km1 = max(k-1, 0);
      for (j = js-1; j <= je+1; j++) {
	u[IDX(i,j,k)] = u[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(im1,j,k)]));
	w[IDX(i,j,k)] = w[IDX(i,j,k)]/(0.5*(d[IDX(i,j,k)] + d[IDX(i,j,km1)]));
	if (fabs(u[IDX(i,j,k)]) > dy[j]/dt) {
	  printf("zeus_y ux warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,im1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		 i,j,k,ie,je,ke,im1);
	  printf("zeus_y ux warning: u,d,d-1=%"GSYM",%"GSYM",%"GSYM"\n", u[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(im1,j,k)]);
	}
	if (fabs(w[IDX(i,j,k)]) > dy[j]/dt) {

	  printf("zeus_y wx warning: i,j,k=%"ISYM",%"ISYM",%"ISYM"  ie,je,ke,km1 = %"ISYM",%"ISYM",%"ISYM",%"ISYM"\n",
		 i,j,k,ie,je,ke,km1);
	  printf("zeus_y wx warning: v,d,d-1=%"GSYM",%"GSYM",%"GSYM"\n", w[IDX(i,j,k)],d[IDX(i,j,k)],d[IDX(i,j,km1)]);
	}
      } // end: loop over j
    } // end: loop over i
  } // end: loop over k

  /* Cleanup */

  for (ic = 0; ic < ncolor; ic++)
    delete [] colstar[ic];

  return SUCCESS;

}
