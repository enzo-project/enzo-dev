/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR ZEUS HYDRO)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     end     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     gridvel - bulk grid velocity (vector of length 3)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (not currently used)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     bottom  - true (1) if this is the lowest level
c     minsupecoef - coefficient for minimum pressure support (0 - not used)
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"


int Zeus_xTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dx[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int NumberOfColours, int colnum[]);

int Zeus_yTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dy[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int NumberOfColours, int colnum[]);

int Zeus_zTransport(float *d, float *e, float *u, float *v, float *w, 
		    int in, int jn, int kn, int rank,
		    int is, int ie, int js, int je, int ks, int ke,
		    float dt, float dz[], float *f1, int bottom,
		    int nsubgrids, long_int GridGlobalStart[],
		    fluxes *SubgridFluxes[], int DensNum, int TENum, 
		    int Vel1Num, int Vel2Num, int Vel3Num, float *BaryonField[],
		    int NumberOfColours, int colnum[]);

int ZeusSource(float *d, float *e, float *u, float *v, float *w, float *p, float *cr, 
	       int in, int jn, int kn, int rank, int igamfield,
	       int is, int ie, int js, int je, int ks, int ke, 
	       float C1, float C2, int ipresfree,
	       float *gamma, float dt, float pmin, float dx[], float dy[], float dz[],
	       int gravity, float *gr_xacc, float *gr_yacc, float *gr_zacc, 
	       int bottom, float minsupecoef, int CRModel, float CRgamma);

int GetUnits (float *DensityUnits, float *LengthUnits,
         float *TemperatureUnits, float *TimeUnits,
         float *VelocityUnits, double *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::ZeusSolver(float *gamma, int igamfield, int nhy, 
		     float dx[], float dy[], float dz[], 
		     int gravity, int NumberOfSubgrids, long_int GridGlobalStart[],
		     fluxes *SubgridFluxes[],
		     int NumberOfColours, int colnum[], int bottom,
		     float minsupecoef)
{

  /*  Locals */

  int i, ie, is, j, je, js, k, ks, ke, n, ixyz, ret;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum;
  float pmin;
  float *d, *e, *u, *v, *w, *cr, *m;

  /* Error Check */

  for (i = 0; i < GridRank; i++)
    if (GridDimension[i] > MAX_ANY_SINGLE_DIRECTION) {
      ENZO_FAIL("ZEUS_MAIN: A grid dimension is too long (increase max_any_single_direction.)\n");
    }

  /* Allocate temporary space for Zeus_Hydro. */

  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float *p = new float[size];
  
  /* Find fields: density, total energy, velocity1-3 and set pointers to them
     Create zero fields for velocity2-3 for low-dimension runs because solver
     assumes they exist. */

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum );
  if (CRModel) {
    if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
      ENZO_FAIL("Cannot Find Cosmic Rays");
    cr = BaryonField[CRNum];
  }
  
  d = BaryonField[DensNum];
  e = BaryonField[TENum];
  u = BaryonField[Vel1Num];
  v = BaryonField[Vel2Num];
  w = BaryonField[Vel3Num];
  if (GridRank < 2) {
    v = new float[size];
    for (i = 0; i < size; i++)
      v[i] = 0;
  }
  if (GridRank < 3) {
    w = new float[size];
    for (i = 0; i < size; i++)
      w[i] = 0;
  }


  /*  Set grid start and end indicies */

  is = GridStartIndex[0];
  js = GridStartIndex[1];
  ks = GridStartIndex[2];
  ie = GridEndIndex[0];
  je = GridEndIndex[1];
  ke = GridEndIndex[2];

  //  If NumberOfGhostZones is set to 4, then use the extra space

  if (is == 4) {
    is = is - 1;
    ie = ie + 1;
  }
  if (js == 4) {
    js = js - 1;
    je = je + 1;
  }
  if (ks == 4) {
    ks = ks - 1;
    ke = ke + 1;
  }

  /* Set minimum pressure (better if it were a parameter) */

  pmin = tiny_number;

  /* Error check */

  for (i = 0; i < size; i++) {
    if (fabs(u[i]) > dx[0]/dtFixed ||
	fabs(v[i]) > dy[0]/dtFixed ||
	fabs(w[i]) > dz[0]/dtFixed) {
      fprintf(stderr, "u,v,w,d,e=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM"  dt=%"GSYM"\n", 
	      u[i],v[i],w[i],d[i],e[i], dx[0], dtFixed);
      ENZO_FAIL("Velocity too fast! (pre-call)\n");
    }
  }
		 
  /*   1) Add source terms */

  if (ZeusSource(d, e, u, v, w, p, cr, 
		 GridDimension[0], GridDimension[1], GridDimension[2],
		 GridRank, igamfield,
		 is, ie, js, je, ks, ke, 
		 ZEUSLinearArtificialViscosity,
		 ZEUSQuadraticArtificialViscosity, PressureFree,
		 gamma, dtFixed, pmin, dx, dy, dz,
		 gravity, AccelerationField[0], AccelerationField[1],
		 AccelerationField[2],
		 bottom, minsupecoef, CRModel, CRgamma) == FAIL) {
    fprintf(stderr, "P(%"ISYM"): Error in ZeusSource on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber,
	    nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
    ENZO_FAIL("Error in ZeusSource!\n");
  }

  /* Error check */
  
  float CRcs = 0.0;
  if (CRmaxSoundSpeed != 0.0){
		  // Get system of units
    float CRsound,DensityUnits,LengthUnits,TemperatureUnits,
          TimeUnits,VelocityUnits,Time;
    double MassUnits;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    CRsound = CRmaxSoundSpeed/VelocityUnits; 
    CRcs = (CRgamma-1.0)/(CRsound*CRsound);
  }

  for (i = 0; i < size; i++) {
    if (fabs(u[i]) > dx[0]/dtFixed ||
	fabs(v[i]) > dy[0]/dtFixed ||
	fabs(w[i]) > dz[0]/dtFixed) {
      fprintf(stderr, "u,v,w,d,e=%"GSYM",%"GSYM",%"GSYM",%"GSYM",%"GSYM"  dx=%"GSYM"  dt=%"GSYM"\n", 
	      u[i],v[i],w[i],d[i],e[i], dx[0], dtFixed);
      ENZO_FAIL("Velocity too fast! (post-call)\n");
    }
  
    /* -- density/TE floor for CR model -- */

    if ( CRModel ){
      if ( CRdensFloor != 0.0 && d[i] < CRdensFloor ) d[i] = CRdensFloor;
      if ( CRcs        != 0.0 && d[i] < CRcs*cr[i]  ) d[i] = CRcs*cr[i];   // Limits sound-speed 
      if ( e[i] < tiny_number*1e-5                  ) e[i] = tiny_number*1e-5;
    } // end cr model if
  } // end i for

  /*  2) Transport step */

  ixyz = nhy % GridRank;
  for(n=ixyz; n <= ixyz+GridRank-1; n++) {

    /* Transport step - x direction */
    
    if ((n % GridRank) == 0)
      ret = Zeus_xTransport(d, e, u, v, w, GridDimension[0], 
			    GridDimension[1], GridDimension[2], GridRank,
			    is, ie, js, je, ks, ke,
			    dtFixed, dx, p, bottom,
			    NumberOfSubgrids, GridGlobalStart,
			    SubgridFluxes, DensNum, TENum,
			    Vel1Num, Vel2Num, Vel3Num, BaryonField,
			    NumberOfColours, colnum);

    /*  Transport step - y direction */

    if ((n % GridRank) == 1 && GridRank > 1)
      ret = Zeus_yTransport(d, e, u, v, w, GridDimension[0], 
			    GridDimension[1], GridDimension[2], GridRank,
			    is, ie, js, je, ks, ke,
			    dtFixed, dy, p, bottom,
			    NumberOfSubgrids, GridGlobalStart,
			    SubgridFluxes, DensNum, TENum,
			    Vel1Num, Vel2Num, Vel3Num, BaryonField,
			    NumberOfColours, colnum);

    /*  Transport step - z direction */

    if ((n % GridRank) == 2 && GridRank > 2)
      ret = Zeus_zTransport(d, e, u, v, w, GridDimension[0], 
			    GridDimension[1], GridDimension[2], GridRank,
			    is, ie, js, je, ks, ke,
			    dtFixed, dz, p, bottom,
			    NumberOfSubgrids, GridGlobalStart,
			    SubgridFluxes, DensNum, TENum,
			    Vel1Num, Vel2Num, Vel3Num, BaryonField,
			    NumberOfColours, colnum);

    if (ret == FAIL) {
      fprintf(stderr, "P(%"ISYM"): Error on ZeusTransport dim=%"ISYM" (Cycle = %"ISYM", dt=%"GSYM")\n", 
	      MyProcessorNumber, n % GridRank, nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
      ENZO_FAIL("Error in ZeusSource!\n");
    }
  
  } // end loop over n

  /* Clean up */

  delete [] p;
  if (GridRank < 2) delete [] v;
  if (GridRank < 3) delete [] w;

  
  
  return SUCCESS;

}
