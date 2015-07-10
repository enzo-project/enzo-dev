/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF STAR PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  April, 2009 by JHW to have multiple types of star 
/              particles
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
#define NO_STAR1
 
#ifdef STAR1
extern "C" void FORTRAN_NAME(star_maker1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *dx, FLOAT *t, float *z, int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff, hydro_method *imethod,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp);
#endif /* STAR1 */
 
extern "C" void FORTRAN_NAME(star_maker2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
		float *mp, float *tdp, float *tcp, float *metalf,
	     int *imetalSNIa, float *metalSNIa, float *metalfSNIa);
 
extern "C" void FORTRAN_NAME(star_maker3mom)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
		 float *mp, float *tdp, float *tcp, float *metalf,
 	     int *imetalSNIa, float *metalSNIa, float *metalfSNIa);

extern "C" void FORTRAN_NAME(star_maker3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *zfield1, float *zfield2,
             float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
		 float *mp, float *tdp, float *tcp, float *metalf,
 	     int *imetalSNIa, float *metalSNIa, float *metalfSNIa);

extern "C" void FORTRAN_NAME(star_maker4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
	         int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf,
 	     int *imetalSNIa, float *metalSNIa, float *metalfSNIa);


 extern "C" void FORTRAN_NAME(star_maker7)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z,
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     		 int *ibuff,
             int *imetal, hydro_method *imethod, float *mintdyn,
             float *odthresh, float *massff, float *smthrest, int *level,
		 int *np, int *npart,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, 
	     FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, 
		   int *type, int *ctype, int *option,
	     int *imetalSNIa, float *metalSNIa, float *metalfSNIa);


extern "C" void FORTRAN_NAME(star_maker5)
  (int *nx, int *ny, int *nz,
   float *d, float *dm, float *temp, float *coolrate, float *u, 
        float *v, float *w, float *cooltime,
   float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
        int *procnum,
   float *d1, float *x1, float *v1, float *t1,
   int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff, 
   int *imetal, hydro_method *imethod, float *mintdyn,
   float *odthresh, float *shdens, 
   float *massff, float *smthrest, int *level, int *np,
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
   float *mp, float *tdp, float *tcp, float *metalf,
   FLOAT *rr_left0, FLOAT *rr_left1, FLOAT *rr_left2, 
   FLOAT *rr_right0, FLOAT *rr_right1, FLOAT *rr_right2,
   int *imetalSNIa, float *metalSNIa, float *metalfSNIa);


int star_maker8(int *nx, int *ny, int *nz, int *size,
		float *d, float *te, float *ge, float *u, float *v, float *w,
		float *bx, float *by, float *bz, float *dt, float *r, float *dx, FLOAT *t, float *z, 
		int *procnum,
		float *d1, float *x1, float *v1, float *t1,
		int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		int *ibuff, 
		int *imethod, int *idual, float *massthresh, int *level, int *np,
		FLOAT *xp, FLOAT *yp, FLOAT *zp, 
		float *up, float *vp, float *wp, float *mp, 
		float *tcp, float *tdp, float *dm, int *type,
		int *npold, FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, 
		float *upold, float *vpold, float *wpold, float *mpold,
		float *tcpold, float *tdpold, float *dmold, 
		float *nx_jet, float *ny_jet, float *nz_jet,
		int *typeold, PINT *idold, int *ctype,
		float *jlrefine, float *temp, float *gamma, float *mu,
		int *nproc, int *nstar);

int star_maker9(int *nx, int *ny, int *nz, int *size,
		float *d, float *te, float *ge, float *u, float *v, float *w,
		float *bx, float *by, float *bz, float *dt, float *r, float *dx, FLOAT *t, float *z, 
		int *procnum,
		float *d1, float *x1, float *v1, float *t1,
		int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		int *ibuff, 
		int *imethod, int *idual, float *massthresh, int *level, int *np,
		FLOAT *xp, FLOAT *yp, FLOAT *zp, 
		float *up, float *vp, float *wp, float *mp, 
		float *tcp, float *tdp, float *dm, int *type,
		int *npold, FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, 
		float *upold, float *vpold, float *wpold, float *mpold,
		float *tcpold, float *tdpold, float *dmold, 
		float *nx_jet, float *ny_jet, float *nz_jet,
		int *typeold, PINT *idold, int *ctype,
		float *jlrefine, float *temp, float *gamma, float *mu,
		int *nproc, int *nstar);

extern "C" void FORTRAN_NAME(star_maker_h2reg)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v, float *w,
                float *cooltime,
             float *dt, float *r, float *metal, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
             int *imetal, hydro_method *imethod, 
	     float *StarFormationEfficiency,
             float *StarFormationNumberDensityThreshold, 
             float *StarFormationMinimumMass, 
             float *MinimumH2FractionForStarFormation,
             int *StochasticStarFormation,
             int *UseSobolevColumn,
             float *SigmaOverR,
             int *AssumeColdWarmPressureBalance,
             float *H2DissociationFlux_MW, 
             float *H2FloorInColdGas,
             float *ColdGasTemperature,
	     int *ran1_init,
             int *level, int *grid_id,
	     int *oldnp, FLOAT *oldxp, FLOAT *oldyp, FLOAT *oldzp,
	     float *oldmp, float *oldtcp, int *oldtype,
	     int *np, FLOAT *xp, FLOAT *yp, FLOAT *zp, 
             float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf);


#ifdef STAR1
extern "C" void FORTRAN_NAME(star_feedback1)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *temp, float *u, float *v,
		       float *w, float *dt, float *r, float *dx,
                       FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
     				 int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *te, float *ge,
		       int *idual);
#endif /* STAR1 */
 
extern "C" void FORTRAN_NAME(star_feedback2)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
	     int *distrad, int *diststep, int *distcells,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn);
 
extern "C" void FORTRAN_NAME(star_feedback3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal, float *zfield1, float *zfield2,
	     int *idual, int *imetal, int *imulti_metals, hydro_method *imethod, 
		       float *dt, float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
	     int *distrad, int *diststep, int *distcells,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn, float *kinf);

extern "C" void FORTRAN_NAME(star_feedback30)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal, float *zfield1, float *zfield2,
	     int *idual, int *imetal, int *imulti_metals, hydro_method *imethod, 
		       float *dt, float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
	     int *distrad, int *diststep, int *distcells,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback5)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn);
 
extern "C" void FORTRAN_NAME(star_feedback4)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback7)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v, 
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
		       int *ibuff, int *level,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int *type,
			float *justburn, int *ctype, float *mbhradius);

extern "C" void FORTRAN_NAME(star_feedback_pn_snia)
  (int *nx, int *ny, int *nz,
   float *d, float *dm, float *te, float *ge, float *u, float *v,
       float *w, float *metal,
   int *idual, int *imetal, hydro_method *imethod, float *dt,
       float *r, float *dx, FLOAT *t, float *z,
   float *d1, float *x1, float *v1, float *t1,
       float *sn_param, float *m_eject, float *yield,
   int *distrad, int *diststep, int *distcells,
       int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff,
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
       float *mp, float *tdp, float *tcp, float *metalf, int *type,
   float *justburn, int *iPN, int *imetalSNIa, float *metalSNIa);

extern "C" void FORTRAN_NAME(pop3_maker)
  (int *nx, int *ny, int *nz, 
   float *d, float *dm, float *h2d, float *temp, 
   float *u, float *v, float *w, 
   float *cooltime, float *dt, float *r, float *metal, float *dx, FLOAT *t, 
   float *z, int *procnum, 
   float *d1, float *x1, float *v1, float *t1, 
   int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
   int *ibuff, int *imetal, hydro_method *imethod, float *h2crit, 
   float *metalcrit, float *odthresh, float *starmass, int *level, int *np, 
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp, 
   float *mp, float *tdp, float *tcp, float *metalf, 
   int *type, int *ctype, float *justburn, int *iradtrans);

extern "C" void PFORTRAN_NAME(pop3_color_maker)
  (int *nx, int *ny, int *nz, 
   float *d, float *dm, float *u, float *v, float *w, 
   float *dt, float *r, float *dx, FLOAT *t, float *z, int *procnum, 
   float *d1, float *x1, float *v1, float *t1, 
   int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff,
   hydro_method *imethod, float *odthresh, int *level, int *np, 
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp, 
   float *mp, float *tcp, int *type, int *ctype);

extern "C" void FORTRAN_NAME(cluster_maker)
  (int *nx, int *ny, int *nz, 
   float *formleft, float *formright,
   float *d, float *dm, float *temp, 
   float *u, float *v, float *w, 
   float *cooltime, float *dt, float *r, float *metal, float *dx, FLOAT *t, 
   float *z, int *procnum, 
   float *d1, float *x1, float *v1, float *t1, 
   int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff, 
   int *imetal, hydro_method *imethod, int *imulti, float *efficiency, 
   float *metalcrit, float *odthresh, float *lifetime, int *level, int *np, 
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp, 
   float *mp, float *tdp, float *tcp, float *metalf, 
   int *type, int *ctype, float *justburn, int *iradtrans,
   int *imetalSNIa, float *metalSNIa, float *metalfSNIa);


int sink_maker(int *nx, int *ny, int *nz, int *size,
             float *d, float *u, float *v, float *w,
             float *dt, float *r, float *dx, FLOAT *t, float *z, 
             int *procnum,
             float *d1, float *x1, float *v1, float *t1,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
     		 int *ibuff, 
	       int *imethod, float *massthresh, int *level, int *np, 
             FLOAT *xp, FLOAT *yp, FLOAT *zp, 
	     float *up, float *vp, float *wp, float *mp, 
		 float *tcp, float *tdp, int *type,
             int *npold, FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, 
	     float *upold, float *vpold, float *wpold, float *mpold,
		 float *tcpold, float *tdpold, int *typeold, int *ctype,
	     float *jlrefine, float *temp, float *JRCT);

int mbh_maker(int *nx, int *ny, int *nz, int *size, float *d, float *u, 
	      float *v, float *w, float *dt, float *r, float *dx, FLOAT *t, 
	      float *d1, float *x1, float *v1, float *t1, 
	      int *nmax, FLOAT *xstart, FLOAT *ystart, 
	      FLOAT *zstart, int *ibuff, 
	      int *level, int *np, FLOAT *xp, FLOAT *yp, 
	      FLOAT *zp, float *up, float *vp, float *wp, float *mp, 
	      float *tcp, float *tdp, int *type, int *ctype);
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
// declaring Geoffrey's Emissivity field prototype
#ifdef EMISSIVITY
  int CalcEmiss(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
	      float *justburn, float *EmissivityArray, float dtLevelAbove);
#endif 
 
int grid::StarParticleHandler(HierarchyEntry* SubgridPointer, int level,
			      float dtLevelAbove, float TopGridTimeStep)
{

  if (!StarParticleCreation && !StarParticleFeedback)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

//  printf("Npart = %d\n", NumberOfParticles);
 
  /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);

  /* initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = DEFAULT_GHOST_ZONES;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  const double m_h = 1.673e-24;

  /* Find Multi-species fields. */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum; 
  if (STARMAKE_METHOD(H2REG_STAR))
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }

  /* If only star cluster formation, check now if we're restricting
     formation in a region. */

  if (StarParticleCreation == (1 << STAR_CLUSTER))
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (StarClusterRegionLeftEdge[dim] >= GridRightEdge[dim] ||
	  StarClusterRegionRightEdge[dim] <= GridLeftEdge[dim])
	return SUCCESS;

  LCAPERF_START("grid_StarParticleHandler");

  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
 
  /* If using MHD, subtract magnetic energy from total energy because 
     density may be modified in star_maker8. */
  
  float *Bfieldx = NULL, *Bfieldy = NULL, *Bfieldz = NULL;
  if (HydroMethod == MHD_RK) {
    Bfieldx = BaryonField[B1Num];
    Bfieldy = BaryonField[B2Num];
    Bfieldz = BaryonField[B3Num];
    for (int n = 0; n < size; n++) {
      float den = BaryonField[DensNum][n];
      float Bx  = BaryonField[B1Num  ][n];
      float By  = BaryonField[B2Num  ][n];
      float Bz  = BaryonField[B3Num  ][n];
      float B2 = Bx*Bx + By*By + Bz*Bz;
      BaryonField[TENum][n] -= 0.5*B2/den;
    }
  }

  if (MultiSpecies > 1) {
      H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
      H2IINum  = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
  }

  /* Find metallicity field and set flag. */
 
  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  /* Set variables to type defines to pass to FORTRAN routines */

  int NormalStarType = PARTICLE_TYPE_STAR;
  int SinkParticleType = PARTICLE_TYPE_MUST_REFINE;
  int SingleStarType = PARTICLE_TYPE_SINGLE_STAR;
  int StarClusterType = PARTICLE_TYPE_CLUSTER;
  int MBHParticleType = PARTICLE_TYPE_MBH;
  int ColorStar = PARTICLE_TYPE_COLOR_STAR;
  int SimpleSource = PARTICLE_TYPE_SIMPLE_SOURCE;

  /* Compute the redshift. */
 
  float zred;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  zred = 1.0*(1.0+InitialRedshift)/a - 1.0;
 
  /* Compute the temperature field. */
 
  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);
 
  /* Get the dark matter field in a usable size for star_maker
     (if level > MaximumGravityRefinementLevel then the dark matter
      field is not valid, so just make it zero - by this time, the
      evolution will be dominated by baryonic matter anyway). */
 
  float *dmfield = new float[size];
  int StartIndex[MAX_DIMENSION], Zero[] = {0,0,0};
  if (level <= MaximumGravityRefinementLevel &&
      GravitatingMassFieldParticles != NULL) {
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      StartIndex[dim] =
      nint((CellLeftEdge[dim][0] - GravitatingMassFieldParticlesLeftEdge[dim])/
	   GravitatingMassFieldParticlesCellSize);
    FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, dmfield,
			 GravitatingMassFieldParticlesDimension,
			 GravitatingMassFieldParticlesDimension+1,
			 GravitatingMassFieldParticlesDimension+2,
			 GridDimension, GridDimension+1, GridDimension+2,
			 Zero, Zero+1, Zero+2,
			 StartIndex, StartIndex+1, StartIndex+2);
  } else
    for (i = 0; i < size; i++)
      dmfield[i] = 0;
 
  /* Convert the species densities into fractional densities (i.e. divide
     by total baryonic density).  At the end we will multiply by the new
     density so that species fractions are maintained. */
 
  for (field = 0; field < NumberOfBaryonFields; field++)
    if ((FieldType[field] >= ElectronDensity && FieldType[field] <= ExtraType1) ||
	FieldType[field] == MetalSNIaDensity)
#ifdef EMISSIVITY
      /* 
         it used to be set to  FieldType[field] < GravPotential if Geoffrey's Emissivity0
         baryons field is used, but no longer needed since it is set to <=ExtraType1
         so the values will scale inside StarParticleHandler 
      */
#endif
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0] +
	    GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	    BaryonField[field][index] /= BaryonField[DensNum][index];
	}

  /* If creating primordial stars, make a total H2 density field */

  float *h2field = NULL;
  if (STARMAKE_METHOD(POP3_STAR)) {
    h2field = new float[size];
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + 
	  GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) 
	  h2field[index] = BaryonField[H2INum][index] + BaryonField[H2IINum][index];
      }
  }

  /* If both metal fields exist, make a total metal field */

  float *MetalPointer;
  float *TotalMetals = NULL;
  int MetallicityField;

  MetallicityField = (MetalNum != -1 || SNColourNum != -1);

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  //printf("Star type \n");
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  /* If we're using physical units for the overdensity threshold and
     it is more than 10 times the cosmic mean, convert it from /cm3 to
     overdensity.  Currently only for Pop II (radiating clusters)/III
     star particles.  */

  float OverDensityThreshold;
  if (PopIIIOverDensityThreshold < 0) {
    OverDensityThreshold = -PopIIIOverDensityThreshold * 1.673e-24 / DensityUnits;
    if (OverDensityThreshold < 1)
      OverDensityThreshold = huge_number;
  }
  else
    OverDensityThreshold = PopIIIOverDensityThreshold;
 
  float CellWidthTemp = float(CellWidth[0][0]);
  float PopIIIMass = (PopIIIInitialMassFunction == TRUE) ? 
    PopIIILowerMassCutoff : PopIIIStarMass;
 
  /* ------------------------------------------------------------------- */
  /* 1) StarParticle creation. */
 
  //  if (StarParticleCreation > 0 && level == MaximumRefinementLevel) {
  if (StarParticleCreation > 0) {
    
    /* Generate a fake grid to keep the particles in. */
 
    grid *tg = new grid;
    tg->GridRank = GridRank;
    tg->ProcessorNumber = ProcessorNumber;
 
    /* Allocate space for new particles. */
 
    int MaximumNumberOfNewParticles = int(0.25*float(size)) + 5;
    tg->AllocateNewParticles(MaximumNumberOfNewParticles);
 
    /* Compute the cooling time. */
 
    float *cooling_time = new float[size];
    this->ComputeCoolingTime(cooling_time);
 
    /* Call FORTRAN routine to do the actual work. */
 
    int NumberOfNewParticlesSoFar = 0;
    int NumberOfNewParticles = 0;
 
#ifdef STAR1
    //    if (StarParticleCreation == 1) {
    if (0) {
      FORTRAN_NAME(star_maker1)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, &HydroMethod,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
          &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0]);
    }
#endif /* STAR1 */
 
    if (STARMAKE_METHOD(NORMAL_STAR)) {

      //---- MODIFIED SF ALGORITHM ("STANDARD VERSION")

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(star_maker2)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], MetalPointer,
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2],
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);

      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    } 

    if (STARMAKE_METHOD(UNIGRID_STAR_MOM)) {

      //---- UNIGRID ALGORITHM (NO JEANS MASS)
      
      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(star_maker3mom)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
       BaryonField[MetalNum+1], BaryonField[MetalNum+2],
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2],
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);

      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    }

    if (STARMAKE_METHOD(UNIGRID_STAR)) {

      //---- UNIGRID ALGORITHM (NO JEANS MASS)
      
      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(star_maker3)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
       BaryonField[MetalNum+1], BaryonField[MetalNum+2],
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles, 
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2],
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);

      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    }

    if (STARMAKE_METHOD(KRAVTSOV_STAR)) {

      //---- KRAVTSOV SF ALGORITHM

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(star_maker4)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], MetalPointer, 
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, 
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime, 
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1], 
          tg->ParticlePosition[2], 
       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
          tg->ParticleVelocity[2], 
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2],
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);


      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    }

    if (STARMAKE_METHOD(POP3_STAR)) {

      //---- POPULATION III (SINGLE STAR)

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(pop3_maker)
	(GridDimension, GridDimension+1, GridDimension+2, BaryonField[DensNum], 
	 dmfield, h2field, temperature, BaryonField[Vel1Num], 
	 BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time, &dtFixed, 
	 BaryonField[NumberOfBaryonFields], MetalPointer, 
	 &CellWidthTemp, &Time, &zred, &MyProcessorNumber, &DensityUnits, 
	 &LengthUnits, &VelocityUnits, &TimeUnits, &MaximumNumberOfNewParticles, 
	 CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
	 &MetallicityField, &HydroMethod, &PopIIIH2CriticalFraction, 
	 &PopIIIMetalCriticalFraction, &OverDensityThreshold, 
	 &PopIIIMass, &level, &NumberOfNewParticles, 
	 tg->ParticlePosition[0], tg->ParticlePosition[1],
	 tg->ParticlePosition[2], tg->ParticleVelocity[0], 
	 tg->ParticleVelocity[1], tg->ParticleVelocity[2], tg->ParticleMass, 
	 tg->ParticleAttribute[1], tg->ParticleAttribute[0], 
	 tg->ParticleAttribute[2], tg->ParticleType, &SingleStarType, 
	 &RadiationData.IntegratedStarFormation, &RadiativeTransfer);

    }

    if (STARMAKE_METHOD(COLORED_POP3_STAR)) {

      //---- POPULATION III (SINGLE STAR)

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      PFORTRAN_NAME(pop3_color_maker)
        (GridDimension, GridDimension+1, GridDimension+2, 
         BaryonField[DensNum], dmfield, 
           BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], 
         &dtFixed, BaryonField[NumberOfBaryonFields], 
           &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
         &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
         &MaximumNumberOfNewParticles, 
           CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
         &HydroMethod, &PopIIIColorDensityThreshold, 
           &level, &NumberOfNewParticles, 
         tg->ParticlePosition[0], tg->ParticlePosition[1], tg->ParticlePosition[2],
         tg->ParticleVelocity[0], tg->ParticleVelocity[1], tg->ParticleVelocity[2],
         tg->ParticleMass, tg->ParticleAttribute[2], tg->ParticleType, &ColorStar);
         
    }

    if (STARMAKE_METHOD(STAR_CLUSTER)) {

      //---- RADIATIVE STELLAR CLUSTERS

      // Convert into a parameter!
      float StarClusterLifeTime = 20e6;  // yr

      // will be assigned lifetime = 4*tdyn
      if (StarClusterUnresolvedModel == TRUE)
	StarClusterLifeTime = FLOAT_UNDEFINED;

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(cluster_maker)
	(GridDimension, GridDimension+1, GridDimension+2, 
	 StarClusterRegionLeftEdge, StarClusterRegionRightEdge,
	 BaryonField[DensNum], dmfield, temperature, 
	 BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], 
	 cooling_time, &dtFixed, BaryonField[NumberOfBaryonFields], 
	 MetalPointer, &CellWidthTemp, 
	 &Time, &zred, &MyProcessorNumber, &DensityUnits, 
	 &LengthUnits, &VelocityUnits, &TimeUnits, &MaximumNumberOfNewParticles, 
	 CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
	 &MetallicityField, &HydroMethod, &MultiSpecies, 
	 &StarClusterFormEfficiency, &PopIIIMetalCriticalFraction, 
	 &OverDensityThreshold, &StarClusterLifeTime, &level, 
	 &NumberOfNewParticles, 
	 tg->ParticlePosition[0], tg->ParticlePosition[1],
	 tg->ParticlePosition[2], tg->ParticleVelocity[0], 
	 tg->ParticleVelocity[1], tg->ParticleVelocity[2], tg->ParticleMass, 
	 tg->ParticleAttribute[1], tg->ParticleAttribute[0], 
	 tg->ParticleAttribute[2], tg->ParticleType, &StarClusterType, 
	 &RadiationData.IntegratedStarFormation, &RadiativeTransfer,
	 &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);


    }

    if (STARMAKE_METHOD(MBH_PARTICLE)) {

      //---- MASSIVE BLACK HOLE PARTICLE 
      //     (particles are put by hand; location picked at MBHInsertLocationFilename, 
      //      once MBH particles are inserted throughout the whole grid hierarchy,
      //      turn off MBH creation --> this is done in EvolveLevel at the bottom of the hierarchy)

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      if (mbh_maker(GridDimension, GridDimension+1, GridDimension+2, &size, 
		    BaryonField[DensNum], BaryonField[Vel1Num],
		    BaryonField[Vel2Num], BaryonField[Vel3Num],
		    &dtFixed, BaryonField[NumberOfBaryonFields],
		    &CellWidthTemp, &Time, 
		    &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
		    &MaximumNumberOfNewParticles, CellLeftEdge[0], 
		    CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
		    &level, &NumberOfNewParticles, tg->ParticlePosition[0], 
		    tg->ParticlePosition[1], tg->ParticlePosition[2], 
		    tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
		    tg->ParticleVelocity[2], tg->ParticleMass, 
		    tg->ParticleAttribute[0], tg->ParticleAttribute[1], 
		    tg->ParticleType, &MBHParticleType) == FAIL) {
	ENZO_FAIL("Error in mbh_maker.");
      }
      
    }

    if (STARMAKE_METHOD(INSTANT_STAR)) {

      //---- MODIFIED SF ALGORITHM (NO-JEANS MASS, NO dt DEPENDENCE, NO STOCHASTIC SF, 
      //                            can reconsider star formation or feedback when MBH exists,
      //                            when in cosmological sim, StarMakerOverDensity is in particles/cc, 
      //                            not the ratio with respect to the DensityUnits, unlike others)

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      // change the unit for StarMakerOverDensity for cosmological run
      if (ComovingCoordinates)
	StarMakerOverDensityThreshold *= m_h / DensityUnits;   

      FORTRAN_NAME(star_maker7)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], MetalPointer,
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones,
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime,
       &StarMakerOverDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles, &NumberOfParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2], 
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleType, &MBHParticleType, &MBHTurnOffStarFormation,
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);

      // make it back to original 
      if (ComovingCoordinates)
	StarMakerOverDensityThreshold /= m_h / DensityUnits;  

      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    } 

    if (STARMAKE_METHOD(SPRINGEL_HERNQUIST_STAR)) {

      //---- Springel & Hernquist 2003 SF algorithm

      float *coolingrate = new float[size];
      for(int coolindex=0; coolindex<size; coolindex++){

        float cgsdensity = BaryonField[DensNum][coolindex]*DensityUnits;
        float *electronguessptr;
        float electronguess = 0.01;
        electronguessptr = &electronguess;
        coolingrate[coolindex] = GadgetCoolingRate
          (log10(temperature[coolindex]), cgsdensity, electronguessptr, zred);
      }

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      FORTRAN_NAME(star_maker5)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, coolingrate,
          BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], MetalPointer,
          &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, CellLeftEdge[0], CellLeftEdge[1],
          CellLeftEdge[2], &GhostZones, 
       &MetallicityField, &HydroMethod, &StarMakerMinimumDynamicalTime, 
       &StarMakerOverDensityThreshold, &StarMakerSHDensityThreshold, &StarMakerMassEfficiency,
       &StarMakerMinimumMass, &level, &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1], 
          tg->ParticlePosition[2], 
       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
          tg->ParticleVelocity[2], 
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
          tg->ParticleAttribute[2],
       &RefineRegionLeftEdge[0], &RefineRegionLeftEdge[1], &RefineRegionLeftEdge[2],
       &RefineRegionRightEdge[0], &RefineRegionRightEdge[1], &RefineRegionRightEdge[2],
       &StarMakerTypeIaSNe, BaryonField[MetalIaNum], tg->ParticleAttribute[3]);


      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;

      delete [] coolingrate;

    }


    /* H2-regulated star formation. */

    if ( STARMAKE_METHOD(H2REG_STAR) && 
	 ( this->MakeStars || !StarFormationOncePerRootGridTimeStep ) ) {

      if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                    HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
      }

      NumberOfNewParticlesSoFar = NumberOfNewParticles;

      /* If StarFormationOncePerRootGridTimeStep, use the root grid
	 timestep for the star particle creation. Stars are only
	 created once per root level time step. */
      float TimeStep;
      if(StarFormationOncePerRootGridTimeStep)
	TimeStep = TopGridTimeStep;
      else
	TimeStep = dtFixed;

      FORTRAN_NAME(star_maker_h2reg)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], BaryonField[HINum], temperature,
       BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
       cooling_time, &TimeStep, BaryonField[NumberOfBaryonFields], 
       MetalPointer, &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &MaximumNumberOfNewParticles, 
       CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
       &MetallicityField, &HydroMethod,
       &H2StarMakerEfficiency,
       &H2StarMakerNumberDensityThreshold, 
       &H2StarMakerMinimumMass,
       &H2StarMakerMinimumH2FractionForStarFormation,
       &H2StarMakerStochastic,
       &H2StarMakerUseSobolevColumn,
       &H2StarMakerSigmaOverR, 
       &H2StarMakerAssumeColdWarmPressureBalance,
       &H2StarMakerH2DissociationFlux_MW, 
       &H2StarMakerH2FloorInColdGas,
       &H2StarMakerColdGasTemperature,
       &ran1_init,
       &level, &ID, 
       &NumberOfParticles,
       ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
       ParticleMass, ParticleAttribute[0], ParticleType, 
       &NumberOfNewParticles,
       tg->ParticlePosition[0], tg->ParticlePosition[1], 
       tg->ParticlePosition[2], 
       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
       tg->ParticleVelocity[2], 
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2]);
      
      for (i = NumberOfNewParticlesSoFar; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;

      /* If StarFormationOncePerRootGridTimeStep, reset
	 this::MakeStars to 0, to prevent further star formation until
	 next top level time step has been taken. */
      if(StarFormationOncePerRootGridTimeStep)
	this->MakeStars = 0;
    }


    /* This creates sink particles which suck up mass off the grid. */

    //    if (STARMAKE_METHOD(SINK_PARTICLE))     printf("   Sink Particle\n"); 
    //if (level == MaximumRefinementLevel)     printf("   Max Refinement\n"); 
    if (STARMAKE_METHOD(SINK_PARTICLE) && level == MaximumRefinementLevel || BigStarFormation > 0) {
      /* Set the density threshold by using the mass in a cell which
	 would have caused another refinement. */
 
      int ihydro = (int) HydroMethod;
      float SinkParticleMassThreshold = huge_number;
      float JeansLengthRefinement = FLOAT_UNDEFINED;
      for (int method = 0; method < MAX_FLAGGING_METHODS; method++) {
	if (CellFlaggingMethod[method] == 2)
	  SinkParticleMassThreshold = MinimumMassForRefinement[method]*
	    pow(RefineBy, level*MinimumMassForRefinementLevelExponent[method]);
	if (CellFlaggingMethod[method] == 6) { 
	  JeansLengthRefinement = RefineByJeansLengthSafetyFactor;
	}
      }

      if(BigStarFormation > 0){
	/* set pointer to the wind direction if wind feedback is used*/

	float *nx_jet = NULL, *ny_jet = NULL, *nz_jet = NULL;
	/*printf("Grid_StarParticleHandler l479 - just made nx_jet etc.\n");*/
	if (StellarWindFeedback) {
	  nx_jet = ParticleAttribute[3];
	  ny_jet = ParticleAttribute[4];
	  nz_jet = ParticleAttribute[5];
	}
	//printf("      test line\n");
	if (star_maker9(GridDimension, GridDimension+1, GridDimension+2, &size, 
			BaryonField[DensNum], BaryonField[TENum], BaryonField[GENum],
			BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
			Bfieldx, Bfieldy, Bfieldz,
			&dtFixed, BaryonField[NumberOfBaryonFields],
			&CellWidthTemp, &Time, &zred, &MyProcessorNumber,
			&DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
			&MaximumNumberOfNewParticles, CellLeftEdge[0], 
			CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
			&ihydro, &DualEnergyFormalism, &SinkParticleMassThreshold, &level, 
			&NumberOfNewParticles, tg->ParticlePosition[0], 
			tg->ParticlePosition[1], tg->ParticlePosition[2], 
			tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
			tg->ParticleVelocity[2], tg->ParticleMass, 
			tg->ParticleAttribute[0], tg->ParticleAttribute[1], tg->ParticleAttribute[2],
			tg->ParticleType, &NumberOfParticles, ParticlePosition[0],
			ParticlePosition[1], ParticlePosition[2], 
			ParticleVelocity[0], ParticleVelocity[1], 
			ParticleVelocity[2], ParticleMass, ParticleAttribute[0], 
			ParticleAttribute[1], ParticleAttribute[2], nx_jet, ny_jet, nz_jet,
			ParticleType, ParticleNumber, &SimpleSource, 
			&JeansLengthRefinement, temperature, &Gamma, &Mu,
			&MyProcessorNumber, &NumberOfStarParticles) == FAIL) {
	  ENZO_FAIL("Error in star_maker9.\n");

	}
      }
      else if(StellarWindFeedback || HydroMethod == MHD_RK || HydroMethod == HD_RK || ProblemType == 107){
	/* set pointer to the wind direction if wind feedback is used*/

	float *nx_jet = NULL, *ny_jet = NULL, *nz_jet = NULL;
	/*printf("Grid_StarParticleHandler l479 - just made nx_jet etc.\n");*/
	if (StellarWindFeedback) {
	  nx_jet = ParticleAttribute[3];
	  ny_jet = ParticleAttribute[4];
	  nz_jet = ParticleAttribute[5];
	}
	//printf("      test line\n");
	if (star_maker8(GridDimension, GridDimension+1, GridDimension+2, &size, 
			BaryonField[DensNum], BaryonField[TENum], BaryonField[GENum],
			BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
			Bfieldx, Bfieldy, Bfieldz,
			&dtFixed, BaryonField[NumberOfBaryonFields],
			&CellWidthTemp, &Time, &zred, &MyProcessorNumber,
			&DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
			&MaximumNumberOfNewParticles, CellLeftEdge[0], 
			CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
			&ihydro, &DualEnergyFormalism, &SinkParticleMassThreshold, &level, 
			&NumberOfNewParticles, tg->ParticlePosition[0], 
			tg->ParticlePosition[1], tg->ParticlePosition[2], 
			tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
			tg->ParticleVelocity[2], tg->ParticleMass, 
			tg->ParticleAttribute[0], tg->ParticleAttribute[1], tg->ParticleAttribute[2],
			tg->ParticleType, &NumberOfParticles, ParticlePosition[0],
			ParticlePosition[1], ParticlePosition[2], 
			ParticleVelocity[0], ParticleVelocity[1], 
			ParticleVelocity[2], ParticleMass, ParticleAttribute[0], 
			ParticleAttribute[1], ParticleAttribute[2], nx_jet, ny_jet, nz_jet,
			ParticleType, ParticleNumber, &SinkParticleType, 
			&JeansLengthRefinement, temperature, &Gamma, &Mu,
			&MyProcessorNumber, &NumberOfStarParticles) == FAIL) {
	  ENZO_FAIL("Error in star_maker8.\n");
	}
      } else {
	  printf("Grid_StarParticleHandler 784 - sink maker called (NOT STAR_MAKER8)\n");
	if (sink_maker(GridDimension, GridDimension+1, GridDimension+2, &size, 
		       BaryonField[DensNum], BaryonField[Vel1Num],
		       BaryonField[Vel2Num], BaryonField[Vel3Num],
		       &dtFixed, BaryonField[NumberOfBaryonFields],
		       &CellWidthTemp, &Time, &zred, &MyProcessorNumber,
		       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
		       &MaximumNumberOfNewParticles, CellLeftEdge[0], 
		       CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
		       &ihydro, &SinkParticleMassThreshold, &level, 
		       &NumberOfNewParticles, tg->ParticlePosition[0], 
		       tg->ParticlePosition[1], tg->ParticlePosition[2], 
		       tg->ParticleVelocity[0], tg->ParticleVelocity[1], 
		       tg->ParticleVelocity[2], tg->ParticleMass, 
		       tg->ParticleAttribute[0], tg->ParticleAttribute[1], 
		       tg->ParticleType, &NumberOfParticles, ParticlePosition[0],
		       ParticlePosition[1], ParticlePosition[2], 
		       ParticleVelocity[0], ParticleVelocity[1], 
		       ParticleVelocity[2], ParticleMass, ParticleAttribute[0], 
		       ParticleAttribute[1], ParticleType, &SinkParticleType, 
		       &JeansLengthRefinement, temperature,
               &JeansRefinementColdTemperature) == FAIL) {
	  ENZO_FAIL("Error in sink_maker.\n");
	}
      }

      /* Delete any merged particles (Mass == FLOAT_UNDEFINED) */
      
      for (int n = 0; n < NumberOfParticles; n++)
	if (ParticleType[n] == SinkParticleType && 
	    ParticleMass[n] == FLOAT_UNDEFINED)
	  NumberOfStarParticles--;

      this->CleanUpMovedParticles();

    } // ENDIF sinks
 
    /* If not set in the above routine, then set the metal fraction to zero. */

    int NoMetallicityAttribute = STARMAKE_METHOD(SINK_PARTICLE);
    if (MetallicityField == FALSE || NoMetallicityAttribute)
      for (i = 0; i < NumberOfNewParticles; i++)
	tg->ParticleAttribute[2][i] = 0.0;
 
    delete [] cooling_time;

      /* Add magnetic energy to total energy with the new density field */
    if (HydroMethod == MHD_RK)
      for (int n = 0; n < size; n++) {
	float den = BaryonField[DensNum][n];
	float Bx  = BaryonField[B1Num  ][n];
	float By  = BaryonField[B2Num  ][n];
	float Bz  = BaryonField[B3Num  ][n];
	float B2 = Bx*Bx + By*By + Bz*Bz;
	BaryonField[TENum][n] += 0.5*B2/den;
      }
 

    /* Move any new particles into their new homes. */
 
    if (NumberOfNewParticles > 0) {
 
      if (debug)
	printf("Grid_StarParticleHandler: New StarParticles = %"ISYM"\n", NumberOfNewParticles);
 
      /* Set the particle numbers.  The correct indices will be assigned in 
	 CommunicationUpdateStarParticleCount in StarParticleFinalize later.*/
 
      for (i = 0; i < NumberOfNewParticles; i++)
 	tg->ParticleNumber[i] = INT_UNDEFINED;
 
      /* Move Particles into this grid (set cell size) using the fake grid. */
 
      tg->NumberOfParticles = NumberOfNewParticles;
      for (dim = 0; dim < GridRank; dim++) {
	tg->CellWidth[dim] = new FLOAT[1];
	tg->CellWidth[dim][0] = CellWidth[dim][0];
      }
      this->MoveAllParticles(1, &tg);
 
    } // end: if (NumberOfNewParticles > 0)

    /* Clean up and keep it quiet. */

    delete tg; // temporary grid

    //    if (debug) printf("StarParticle: end\n");
 
  }
#ifdef EMISSIVITY
    if (StarMakerEmissivityField > 0) {

      /* where values of the emissivity field is calculated */

      int EtaNum = FindField(Emissivity0, FieldType, NumberOfBaryonFields);

      CalcEmiss(GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
	  ParticleAttribute[2], &RadiationData.IntegratedStarFormation, 
       BaryonField[EtaNum], dtLevelAbove);
    }
#endif 

  /* ------------------------------------------------------------------- */
  /* 2) StarParticle feedback. */
 
#ifdef STAR1
  //if (StarParticleFeedback == 1) {
  if (0) {

    //---- THIS IS THE ORIGINAL ENZO STAR FORMATION ALG.

      FORTRAN_NAME(star_feedback1)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, temperature,
          BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num],
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0], 
        BaryonField[TENum], BaryonField[GENum], &DualEnergyFormalism);
 
  } // end: if (StarParticleFeedback == 1)
#endif /* STAR1 */
 
  if (STARFEED_METHOD(NORMAL_STAR)) {

    //---- THIS IS THE MODIFIED STAR FORMATION ALGORITHM
 
      FORTRAN_NAME(star_feedback2)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], MetalPointer,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield, &StarFeedbackDistRadius, &StarFeedbackDistCellStep, 
       &StarFeedbackDistTotalCells,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
       ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation);
 
  } // end: if NORMAL_STAR
 
  if (STARFEED_METHOD(UNIGRID_STAR)) {

    //---- UNIGRID (NON-JEANS MASS) VERSION
 
      FORTRAN_NAME(star_feedback3)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
          BaryonField[MetalNum+1], BaryonField[MetalNum+2],
       &DualEnergyFormalism, &MetallicityField, &MultiMetals, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield, &StarFeedbackDistRadius, &StarFeedbackDistCellStep, 
       &StarFeedbackDistTotalCells,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
       ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation,
       &StarFeedbackKineticFraction);
 
  } // end: if UNIGRID_STAR

  if (STARFEED_METHOD(UNIGRID_STAR_ORIG)) {

    //---- UNIGRID (NON-JEANS MASS) VERSION
 
      FORTRAN_NAME(star_feedback30)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
          BaryonField[MetalNum+1], BaryonField[MetalNum+2],
       &DualEnergyFormalism, &MetallicityField, &MultiMetals, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield, &StarFeedbackDistRadius, &StarFeedbackDistCellStep, 
       &StarFeedbackDistTotalCells,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
       ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation);
 
  } // end: if UNIGRID_STAR
 
  if (STARFEED_METHOD(KRAVTSOV_STAR)) {  

    //---- KRAVTSOV STAR FORMATION ALGORITHM

      FORTRAN_NAME(star_feedback4)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, 
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], MetalPointer,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod, 
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp, 
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction, 
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1], 
          ParticlePosition[2], 
       ParticleVelocity[0], ParticleVelocity[1], 
          ParticleVelocity[2], 
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation);

  } // end: if KRAVSTOV STAR

  if (STARFEED_METHOD(INSTANT_STAR)) {

    // check whether the particles are correctly located in grids 
    // Ji-hoon Kim in Nov.2009

#define NO_PARTICLE_IN_GRID_CHECK   

#ifdef PARTICLE_IN_GRID_CHECK
    int xindex, yindex, zindex;
    for (i = 0; i < NumberOfParticles; i++) {
      
      xindex = (int)((ParticlePosition[0][i] - CellLeftEdge[0][0]) / CellWidthTemp);
      yindex = (int)((ParticlePosition[1][i] - CellLeftEdge[1][0]) / CellWidthTemp); 
      zindex = (int)((ParticlePosition[2][i] - CellLeftEdge[2][0]) / CellWidthTemp); 

      if (xindex < 0 || xindex > GridDimension[0] || 
	  yindex < 0 || yindex > GridDimension[1] || 
	  zindex < 0 || zindex > GridDimension[2])
	fprintf(stdout, "particle out of grid (C level); xind, yind, zind, level = %d, %d, %d, %d\n",
		xindex, yindex, zindex, level); 
    }
#endif

    //---- MODIFIED SF ALGORITHM (NO-JEANS MASS, NO dt DEPENDENCE, NO STOCHASTIC SF)

      double pc = 3.086e18;
      float mbhradius = MBHFeedbackThermalRadius * pc / LengthUnits; 
 
      FORTRAN_NAME(star_feedback7)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], MetalPointer,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield,
       &NumberOfParticles,
       CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, &level,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2], 
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation, 
       &MBHParticleType, &mbhradius);
 
  } 


  if (STARFEED_METHOD(SPRINGEL_HERNQUIST_STAR)) {  

    //---- SPRINGEL & HERNQUIST ALGORITHM

      FORTRAN_NAME(star_feedback5)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield, 
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], MetalPointer,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod, 
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp, 
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction, 
          &StarMetalYield,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1], 
          ParticlePosition[2], 
       ParticleVelocity[0], ParticleVelocity[1], 
          ParticleVelocity[2], 
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
          ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation);
      //fprintf(stderr,"After feedback is called");

  } // end: if (StarParticleFeedback == 5)

  if (StarMakerTypeIaSNe == 1 || StarMakerPlanetaryNebulae == 1) {

      FORTRAN_NAME(star_feedback_pn_snia)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], MetalPointer,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, BaryonField[NumberOfBaryonFields], &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          &StarEnergyToThermalFeedback, &StarMassEjectionFraction,
          &StarMetalYield, &StarFeedbackDistRadius, &StarFeedbackDistCellStep, 
       &StarFeedbackDistTotalCells,
       &NumberOfParticles,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
       ParticleAttribute[2], ParticleType, &RadiationData.IntegratedStarFormation,
       &StarMakerPlanetaryNebulae, &StarMakerTypeIaSNe, BaryonField[MetalIaNum]);

  }

  /* Convert the species back from fractional densities to real densities. */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if ((FieldType[field] >= ElectronDensity && FieldType[field] <= ExtraType1) ||
	FieldType[field] == MetalSNIaDensity) {
#ifdef EMISSIVITY
      /* 
         it used to be set to  FieldType[field] < GravPotential if Geoffrey's Emissivity0
         baryons field is used, but no longer needed since it is set to <=ExtraType1
         so the values will scale inside StarParticleHandler 
      */
#endif
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0] +
	    GridStartIndex[0];
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	    BaryonField[field][index] *= BaryonField[DensNum][index];
	  }
	}
      }
    }
  }
 
  /* Clean up. */
 
  delete [] h2field;
  delete [] TotalMetals;
  delete [] temperature;
  delete [] dmfield;
  delete [] BaryonField[NumberOfBaryonFields];   // refinement flag field
  BaryonField[NumberOfBaryonFields] = NULL;
 
  //if (debug) printf("StarParticle: end\n");


  LCAPERF_STOP("grid_StarParticleHandler");
  return SUCCESS;
}
