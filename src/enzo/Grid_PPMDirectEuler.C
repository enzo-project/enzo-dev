/***********************************************************************
/
/  GRID CLASS (Perform Direct Euler integration)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  converted from F77 to C++ by Ian McGreer
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// PGI doesn't like Msafeptr, so if it gets turned on, we turn it off here
#pragma global nosafeptr=all

// Wrapper for direct euler hydro solve
//

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(feuler_sweep)
  (int *iter, float *v, float *w,
   int *in, int *jn, int *kn, 
   int *d0n, int *d1n, int *d2n, int *d2nactive,
   int *grav, int *idual, float *eta1, float *eta2,
   int *d0start, int *d0end, int *d1start, int *d1end,
   float *gamma, float *pmin, float *dt, 
   float dx[], float dy[], float dz[],
   int *diff, int *flatten, int *steepen, int *ipresfree, 
   int *ncolour, int *sweep_dir,
   float *dslice, float *eslice, 
   float *uslice, float *vslice, float *wslice,
   float *grslice, float *geslice, float *colslice,
   float *df, float *ef, float *uf, float *vf, float *wf,
   float *gef, float *colf);



int grid::PPMDirectEuler
( int       CycleNumber, 
  int       NumberOfSubgrids, 
  fluxes  * SubgridFluxes[], 
  float   * CellWidthTemp[], 
  long_int  GridGlobalStart[], 
  int       GravityOn,
  int       NumberOfColours, 
  int       colnum[] )
{
  int n, dim, ldim, iter, ixyz;

  /* allocate temporary space for PPM solver (enough to fit 31 of the
       largest possible 2d slices plus 4*NumberOfColours). */

  int tempsize = max(max(GridDimension[0]*GridDimension[1], 
			 GridDimension[1]*GridDimension[2]),
                         GridDimension[2]*GridDimension[0]);

  float *temp = new float[tempsize*(31+NumberOfColours*4)];

  //  Loop over directions, using a Strang-type splitting

  ixyz = CycleNumber % GridRank;

  for (n=ixyz; n<=ixyz+GridRank-1; n++) {

    dim = n % GridRank; // find which direction to do for this sweep

    if (GridEndIndex[dim] - GridStartIndex[dim] > 0) {

      ldim = (dim+2)%3; // compute dimension perpindicular to the slice

      /* Loop over two-dimensional slices. */

//      if (debug) {
//	printf("n=%"ISYM"  dim=%"ISYM"  ldim=%"ISYM"  CycleNumber=%"ISYM"\n", 
//	       n, dim, ldim, CycleNumber);
//      }

      for (iter=0; iter < GridDimension[ldim]; iter++) {
        this->euler_sweep(dim, iter, CycleNumber, NumberOfSubgrids,
	                  SubgridFluxes, CellWidthTemp, 
			  GridGlobalStart, GravityOn,
	                  NumberOfColours, colnum, temp, tempsize);
      }
    }
  }

  /* Clean up. */

  delete [] temp;

  return SUCCESS;
}

//========================================================================
// MACRO DEFINITIONS
//========================================================================

// ************************************************** 
//
// NOTE: "NO_ALIASES" should be defined if compiler optimizations
// generate incorrect code.  In particular, the PGI compilers on
// kraken.nics.utk.edu with -O2 or higher optimization.  The main
// difference is in how indices for the main problem are computed in
// loops in the "GET_INDEX_I" macro.
//
// **************************************************

#ifdef NO_ALIASES

#   define INDEXL(ixl,iyl,izl,nx,ny,nz) \
  (dim==0 ?  ((ixl) + (nx)*((iyl) + (ny)*(izl))) : \
  (dim==1 ?  ((izl) + (nx)*((ixl) + (ny)*(iyl))) : \
             ((iyl) + (nx)*((izl) + (ny)*(ixl)))))

#   define GET_INDEX_I int i = INDEXL(ixl,iyl,izl,nx,ny,nz);

#else

#   define INDEX(ix,iy,iz,nx,ny,nz) ((ix) + (nx)*((iy) + (ny)*(iz)))

#   define GET_INDEX_I int i = INDEX(ix,iy,iz,nx,ny,nz);

#endif

#define DENSNUM 0
#define   TENUM 1
#define VEL1NUM 2
#define VEL2NUM 3
#define VEL3NUM 4
#define  ACCNUM 5
#define   GENUM 6
#define  COLNUM 7
#define    NVAR 8

// COMMENT ME!

#define ATMP(b,v)           &temp[TMPINDEX(b,v)]

// COMMENT ME!

#define SLICE(var,d0,d1)    TMP(0,var,nxl,d0,d1)

// COMMENT ME!

#define TMP(b,v,nd0,d0,d1)   temp[TMPINDEX(b,v)+(nd0)*(d1)+(d0)]

// array offset to divide slice memory storage from field memory storage.
// the total number of slices = #variables + #colors

#define TMPINDEX(b,v)       (((b)*(NVAR+NumberOfColours)+v)*tempsize)

// COMMENT ME!

#define COLOR(ci,d0,d1)      temp[TMPINDEX(0,COLNUM) + ((ci)*nyl+(d1))*nxl+(d0)]

// COMMENT ME!

#define LFACE(n) \
  (SubgridFluxes[n]->LeftFluxStartGlobalIndex[perm[0]][perm[0]]  \
    - GridGlobalStart[perm[0]])

// COMMENT ME!

#define RFACE(n)                                         \
  (SubgridFluxes[n]->RightFluxStartGlobalIndex[perm[0]][perm[0]] \
    - GridGlobalStart[perm[0]] + 1)

// COMMENT ME!

#define LEFTFLUX(var,n,off)  SubgridFluxes[n]->LeftFluxes[var][perm[0]][off]

// COMMENT ME!

#define RIGHTFLUX(var,n,off) SubgridFluxes[n]->RightFluxes[var][perm[0]][off]

// COMMENT ME!

#define LEFTFLUXTMP(var,n,off)   TMP(1,var,nxl,LFACE(n),off)

// COMMENT ME!

#define RIGHTFLUXTMP(var,n,off)  TMP(1,var,nxl,RFACE(n),off)

// COMMENT ME!

#define LCOLFLUXTMP(ci,n,off)  temp[TMPINDEX(1,COLNUM) + ((ci)*nyl+(off))*nxl+LFACE(n)]

// COMMENT ME!

#define RCOLFLUXTMP(ci,n,off)  temp[TMPINDEX(1,COLNUM) + ((ci)*nyl+(off))*nxl+RFACE(n)]


/***********************************************************************
/
/  euler_sweep: Perform one directional sweep of the PPMDE hydro method
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  converted to C++ by Ian McGreer
/  modified2:  refactored by James Bordner to remove unnecessary messy pointers,
/              clarify code intent, and optionally bypass variable aliasing 
/              that confuses some compiler optimizations
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

int grid::euler_sweep
( int       dim, 
  int       iter, 
  int       CycleNumber, 
  int       NumberOfSubgrids, 
  fluxes  * SubgridFluxes[], 
  float   * CellWidthTemp[],
  long_int  GridGlobalStart[], 
  int       GravityOn,
  int       NumberOfColours, 
  int       colnum[], 
  float   * temp, 
  int       tempsize )
{
  int perm[3];
  int o1, o2, ic, d1start, d2start, d1end, d2end, nactive;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int vel1, vel2, vel3;
  int subgrid;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  // Vel1 is the in-sweep velocity dimension, while vel2 and vel3 are
  // the perpindicular velocities.

  switch (dim) {
  case 0: vel1 = Vel1Num;  vel2 = Vel2Num;  vel3 = Vel3Num;  break;
  case 1: vel1 = Vel2Num;  vel2 = Vel3Num;  vel3 = Vel1Num;  break;
  case 2: vel1 = Vel3Num;  vel2 = Vel1Num;  vel3 = Vel2Num;  break;
  }

  // Initialize permutation array perm[] that maps main problem axes
  // to local work axes

  perm[0] = (dim + 0) % 3;
  perm[1] = (dim + 1) % 3;
  perm[2] = (dim + 2) % 3;


#ifdef NO_ALIASES

  // Initialize LOCAL WORK grid indices

  int ixl;
  int iyl;
  int izl;

  // Initialize MAIN PROBLEM grid indices

  /* ix,iy,iz not required */

#else

  // Declare permutation array i3[] for aliasing local work indices
  // with main problem indices.  Note that i3 itself is never directly
  // accessed--it ONLY serves to alias (ix,iy,iz) with some permutation
  // of (ixl,iyl,izl).

  int i3[3];

  // Initialize LOCAL WORK GRID indices

  int & ixl = i3[perm[0]]; // *** NOTE ALIAS ***
  int & iyl = i3[perm[1]]; // *** NOTE ALIAS ***
  int & izl = i3[perm[2]]; // *** NOTE ALIAS ***

  // Initialize MAIN PROBLEM GRID indices

  int & ix  = i3[0];       // *** NOTE ALIAS ***
  int & iy  = i3[1];       // *** NOTE ALIAS ***
  int & iz  = i3[2];       // *** NOTE ALIAS ***

#endif

  // Initialize LOCAL WORK GRID extents

  int nxl = GridDimension[perm[0]];
  int nyl = GridDimension[perm[1]];
  int nzl = GridDimension[perm[2]];
  
  // Initialize MAIN PROBLEM GRID extents

  int nx = GridDimension[0];
  int ny = GridDimension[1];
  int nz = GridDimension[2];

  // Initialize the slice index

  izl = iter;

  /* Copy a slice worth of data from field to slice */

  for (iyl=0; iyl < nyl; iyl++) {
    for (ixl=0; ixl < nxl; ixl++) {
      GET_INDEX_I;
      SLICE(DENSNUM,ixl,iyl) = BaryonField[DensNum][i];
      SLICE(  TENUM,ixl,iyl) = BaryonField[  TENum][i];
      SLICE(VEL1NUM,ixl,iyl) = BaryonField[   vel1][i];
    }

    /* copy 2-velocity slice (or set equal to zero if not
       present). (note that we need to be careful to set it to zero if
       the off-sweep dimension is larger than our rank). */

    if (GridRank > perm[1]) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	SLICE(VEL2NUM,ixl,iyl) = BaryonField[   vel2][i];
      }
    } else {
      for (ixl=0; ixl < nxl; ixl++) {
	SLICE(VEL2NUM,ixl,iyl) = 0;
      }
    }

    /* copy 3-velocity slice (or set equal to zero if not present). */

    if (GridRank > perm[2]) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	SLICE(VEL3NUM,ixl,iyl) = BaryonField[   vel3][i];
      }
    } else {
      for (ixl=0; ixl < nxl; ixl++) {
	SLICE(VEL3NUM,ixl,iyl) = 0;
      }
    }

    /* Copy appropriate acceleration field. */

    if (GravityOn) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
        SLICE(ACCNUM,ixl,iyl) = AccelerationField[dim][i];
      }
    }

    /* Copy Gas energy field if necessary. */

    if (DualEnergyFormalism) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
        SLICE(GENUM,ixl,iyl) = BaryonField[GENum][i];
      }
    }

    /* Copy colors, if needed. */

    for (ic = 0; ic < NumberOfColours; ic++) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	COLOR(ic,ixl,iyl) = BaryonField[colnum[ic]][i];
      }
    }

  } // end: loop over iyl

  // number of active slices:

  nactive = GridEndIndex[perm[2]] - GridStartIndex[perm[2]] + 1; 

  /* Set start and end index on off-sweep direction (do all). */

  int OffSweepStartIndex = 0;
  int OffSweepEndIndex   = nyl - 1;
  int SweepDirection     = perm[0] + 1; // add one to match definition in calcdiss.src:


  /*  Set minimum pressure (better if it were a parameter) */

  float pmin = 1e-20;

  /* Call fortran routine to do actual work. */

  FORTRAN_NAME(feuler_sweep)
    (&iter, BaryonField[vel2], BaryonField[vel3],
     &nx, &ny, &nz,
     &nxl, &nyl, &nzl, &nactive,
     &GravityOn, &DualEnergyFormalism, 
     &DualEnergyFormalismEta1, &DualEnergyFormalismEta2,
     &GridStartIndex[perm[0]], &GridEndIndex[perm[0]], 
     &OffSweepStartIndex, &OffSweepEndIndex,
     &Gamma, &pmin, &dtFixed, 
     CellWidthTemp[perm[0]], CellWidthTemp[perm[1]], CellWidthTemp[perm[2]],
     &PPMDiffusionParameter, &PPMFlatteningParameter, 
     &PPMSteepeningParameter, &PressureFree, &NumberOfColours, &SweepDirection,
     ATMP(0,DENSNUM), ATMP(0,TENUM),
     ATMP(0,VEL1NUM), ATMP(0,VEL2NUM), ATMP(0,VEL3NUM), 
     ATMP(0,ACCNUM),  ATMP(0,GENUM),   ATMP(0,COLNUM),
     ATMP(1,DENSNUM), ATMP(1,TENUM),
     ATMP(1,VEL1NUM), ATMP(1,VEL2NUM), ATMP(1,VEL3NUM), 
     ATMP(1,GENUM),   ATMP(1,COLNUM)
     );

  /* Now copy fluxes computed in this slice into the correct flux storage
     around the subgrids.
     To do this, check this slice against the list of subgrids 
     (all subgrid quantities are zero based). */

  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {

    /* Compute the start and end indicies of the subgrid. */
    
    d1start = SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[perm[0]][perm[1]]
      - GridGlobalStart[perm[1]];
    d1end = SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[perm[0]][perm[1]]
      - GridGlobalStart[perm[1]];
    d2start = SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[perm[0]][perm[2]]
      - GridGlobalStart[perm[2]];
    d2end = SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[perm[0]][perm[2]] 
      - GridGlobalStart[perm[2]];

    /* If this slice intersects the subgrid, then loop over the data. */

    if (izl >= d2start && izl <= d2end) {
      for (iyl = d1start; iyl <= d1end; iyl++) {
	o1 = iyl;
	if (dim == 1) {
	  o2 = (izl-d2start) + (iyl-d1start)*(d2end-d2start+1);
	} else {
	  o2 = (iyl-d1start) + (izl-d2start)*(d1end-d1start+1);
	}

	/* Copy the left and right fluxes from the temporary storage to
	   the appropriate subgrid face. */

        LEFTFLUX( DensNum, subgrid, o2) = LEFTFLUXTMP( DENSNUM, subgrid, o1);
        RIGHTFLUX(DensNum, subgrid, o2) = RIGHTFLUXTMP(DENSNUM, subgrid, o1);
        LEFTFLUX(   TENum, subgrid, o2) = LEFTFLUXTMP(   TENUM, subgrid, o1);
        RIGHTFLUX(  TENum, subgrid, o2) = RIGHTFLUXTMP(  TENUM, subgrid, o1);
        LEFTFLUX(    vel1, subgrid, o2) = LEFTFLUXTMP( VEL1NUM, subgrid, o1);
        RIGHTFLUX(   vel1, subgrid, o2) = RIGHTFLUXTMP(VEL1NUM, subgrid, o1);
        if (GridEndIndex[1] - GridStartIndex[1] > 0) {
          LEFTFLUX(    vel2, subgrid, o2) = LEFTFLUXTMP( VEL2NUM, subgrid, o1);
          RIGHTFLUX(   vel2, subgrid, o2) = RIGHTFLUXTMP(VEL2NUM, subgrid, o1);
	}
        if (GridEndIndex[2] - GridStartIndex[2] > 0) {
          LEFTFLUX(    vel3, subgrid, o2) = LEFTFLUXTMP( VEL3NUM, subgrid, o1);
          RIGHTFLUX(   vel3, subgrid, o2) = RIGHTFLUXTMP(VEL3NUM, subgrid, o1);
	}
	if (DualEnergyFormalism) {
          LEFTFLUX( GENum, subgrid, o2) = LEFTFLUXTMP( GENUM, subgrid, o1);
          RIGHTFLUX(GENum, subgrid, o2) = RIGHTFLUXTMP(GENUM, subgrid, o1);
	}
        for (ic=0; ic<NumberOfColours; ic++) {
          LEFTFLUX( colnum[ic], subgrid, o2) = LCOLFLUXTMP(ic, subgrid, o1);
          RIGHTFLUX(colnum[ic], subgrid, o2) = RCOLFLUXTMP(ic, subgrid, o1);
        }
      }
    }
  }

  /* Copy data back from slice to field. */

  for (iyl=0; iyl < nyl; iyl++) {
    for (ixl=0; ixl < nxl; ixl++) {
      GET_INDEX_I;
      BaryonField[DensNum][i] = SLICE(DENSNUM,ixl,iyl);
      BaryonField[  TENum][i] = SLICE(  TENUM,ixl,iyl);
      BaryonField[   vel1][i] = SLICE(VEL1NUM,ixl,iyl);
    }

    if (GridRank > perm[1]) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	BaryonField[   vel2][i] = SLICE(VEL2NUM,ixl,iyl);
      }
    }

    if (GridRank > perm[2]) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	BaryonField[   vel3][i] = SLICE(VEL3NUM,ixl,iyl);
      }
    }

    if (DualEnergyFormalism) {

      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	BaryonField[GENum][i] = SLICE(GENUM,ixl,iyl);
      }
    }

    for (ic = 0; ic < NumberOfColours; ic++) {
      for (ixl=0; ixl < nxl; ixl++) {
	GET_INDEX_I;
	BaryonField[colnum[ic]][i] = COLOR(ic,ixl,iyl);
      }
    }

  } // end: loop over iyl

  return SUCCESS;
}
