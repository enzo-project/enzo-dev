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
#include "StarParticleData.h"

#define  PROTONMASS  1.6726e-24

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
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
             float *mp, float *tdp, float *tcp, float *metalf);
 
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
             float *mp, float *tdp, float *tcp, float *metalf);
 
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
             float *mp, float *tdp, float *tcp, float *metalf);

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
                 int *np,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
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
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
			float *justburn);
 
extern "C" void FORTRAN_NAME(star_feedback3)(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal, float *zfield1, float *zfield2,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
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
             float *mp, float *tdp, float *tcp, float *metalf, 
			float *justburn);

extern "C" void FORTRAN_NAME(star_feedback7)(int *nx, int *ny, int *nz,
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
			float *justburn);

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

extern "C" void FORTRAN_NAME(pop3_color_maker)
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
   int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, 
   int *ibuff, int *imetal, hydro_method *imethod, float *efficiency, 
   float *metalcrit, float *odthresh, float *lifetime, int *level, int *np, 
   FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp, 
   float *mp, float *tdp, float *tcp, float *metalf, 
   int *type, int *ctype, float *justburn, int *iradtrans);

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
	     float *jlrefine, float *temp);
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
 
 
int grid::StarParticleHandler(HierarchyEntry* SubgridPointer, int level)
{

  if (!StarParticleCreation && !StarParticleFeedback)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* First, set under_subgrid field */

  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);

  /* initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = DEFAULT_GHOST_ZONES;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,H2INum, H2IINum;
 
  /* If only star cluster formation, check now if we're restricting
     formation in a region. */

  if (StarParticleCreation == (1 << STAR_CLUSTER))
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (StarClusterRegionLeftEdge[dim] >= GridRightEdge[dim] ||
	  StarClusterRegionRightEdge[dim] <= GridLeftEdge[dim])
	return SUCCESS;

  JBPERF_START("grid_StarParticleHandler");

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
     density may be modified in star_maker3. */
  
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
 
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  int UseColour = FALSE, SNColourNum;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) 
      != -1)
    UseColour = TRUE;
  else
    SNColourNum = 0;

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = max(MetallicityField, UseColour);

  /* Set variables to type defines to pass to FORTRAN routines */

  int NormalStarType = PARTICLE_TYPE_STAR;
  int SingleStarType = PARTICLE_TYPE_SINGLE_STAR;
  int StarClusterType = PARTICLE_TYPE_CLUSTER;
  int SinkParticleType = PARTICLE_TYPE_MUST_REFINE;
  int ColorStar = PARTICLE_TYPE_COLOR_STAR;

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
    if (FieldType[field] >= ElectronDensity && FieldType[field] <= ExtraType1 )
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
 
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
 
  float CellWidthTemp = float(CellWidth[0][0]);
 
  /* ------------------------------------------------------------------- */
  /* 1) StarParticle creation. */
 
  //  if (StarParticleCreation > 0 && level == MaximumRefinementLevel) {
  if (StarParticleCreation > 0) {
 
    /* Generate a fake grid to keep the particles in. */
 
    //    if (debug) printf("StarParticle: start\n");
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
 
    int NumberOfNewParticles = 0;
 
    if (debug && NumberOfNewParticles > 0) {
       fprintf(stderr, "StarParticle: New StarParticles = "
	       "%"ISYM"\n", NumberOfNewParticles);
    }

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

      FORTRAN_NAME(star_maker2)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
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
          tg->ParticleAttribute[2]);

      for (i = 0; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    } 

    if (STARMAKE_METHOD(UNIGRID_STAR)) {

      //---- UNIGRID ALGORITHM (NO JEANS MASS)

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
          tg->ParticleAttribute[2]);

      for (i = 0; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    }

    if (STARMAKE_METHOD(KRAVTSOV_STAR)) {

      //---- KRAVTSOV SF ALGORITHM

      FORTRAN_NAME(star_maker4)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum], 
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
          tg->ParticleAttribute[2]);

      for (i = 0; i < NumberOfNewParticles; i++)
          tg->ParticleType[i] = NormalStarType;
    }

    if (STARMAKE_METHOD(POP3_STAR)) {

      //---- POPULATION III (SINGLE STAR)

      FORTRAN_NAME(pop3_maker)
	(GridDimension, GridDimension+1, GridDimension+2, BaryonField[DensNum], 
	 dmfield, h2field, temperature, BaryonField[Vel1Num], 
	 BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time, &dtFixed, 
	 BaryonField[NumberOfBaryonFields], BaryonField[MetalNum], 
	 &CellWidthTemp, &Time, &zred, &MyProcessorNumber, &DensityUnits, 
	 &LengthUnits, &VelocityUnits, &TimeUnits, &MaximumNumberOfNewParticles, 
	 CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
	 &MetallicityField, &HydroMethod, &PopIIIH2CriticalFraction, 
	 &PopIIIMetalCriticalFraction, &PopIIIOverDensityThreshold, 
	 &PopIIIStarMass, &level, &NumberOfNewParticles, 
	 tg->ParticlePosition[0], tg->ParticlePosition[1],
	 tg->ParticlePosition[2], tg->ParticleVelocity[0], 
	 tg->ParticleVelocity[1], tg->ParticleVelocity[2], tg->ParticleMass, 
	 tg->ParticleAttribute[1], tg->ParticleAttribute[0], 
	 tg->ParticleAttribute[2], tg->ParticleType, &SingleStarType, 
	 &RadiationData.IntegratedStarFormation, &RadiativeTransfer);

    }

    if (STARMAKE_METHOD(COLORED_POP3_STAR)) {
      FORTRAN_NAME(pop3_color_maker)
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
      float StarClusterLifeTime = 10e6;  // yr

      FORTRAN_NAME(cluster_maker)
	(GridDimension, GridDimension+1, GridDimension+2, 
	 StarClusterRegionLeftEdge, StarClusterRegionRightEdge,
	 BaryonField[DensNum], dmfield, temperature, 
	 BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], 
	 cooling_time, &dtFixed, BaryonField[NumberOfBaryonFields], 
	 BaryonField[MetalNum], &CellWidthTemp, 
	 &Time, &zred, &MyProcessorNumber, &DensityUnits, 
	 &LengthUnits, &VelocityUnits, &TimeUnits, &MaximumNumberOfNewParticles, 
	 CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones, 
	 &MetallicityField, &HydroMethod, &StarClusterFormEfficiency,
	 &PopIIIMetalCriticalFraction, &PopIIIOverDensityThreshold, 
	 &StarClusterLifeTime, &level, &NumberOfNewParticles, 
	 tg->ParticlePosition[0], tg->ParticlePosition[1],
	 tg->ParticlePosition[2], tg->ParticleVelocity[0], 
	 tg->ParticleVelocity[1], tg->ParticleVelocity[2], tg->ParticleMass, 
	 tg->ParticleAttribute[1], tg->ParticleAttribute[0], 
	 tg->ParticleAttribute[2], tg->ParticleType, &StarClusterType, 
	 &RadiationData.IntegratedStarFormation, &RadiativeTransfer);

    }

    /* This creates sink particles which suck up mass off the grid. */

    if (STARMAKE_METHOD(SINK_PARTICLE) && level == MaximumRefinementLevel) {

      /* Set the density threshold by using the mass in a cell which
	 would have caused another refinement. */

      int ihydro = (int) HydroMethod;
      float SinkParticleMassThreshold = huge_number;
      float JeansLengthRefinement = FLOAT_UNDEFINED;
      for (int method = 0; method < MAX_FLAGGING_METHODS; method++) {
	if (CellFlaggingMethod[method] == 2)
	  SinkParticleMassThreshold = MinimumMassForRefinement[method]*
	    pow(RefineBy, level*MinimumMassForRefinementLevelExponent[method]);
	if (CellFlaggingMethod[method] == 6)
	  JeansLengthRefinement = RefineByJeansLengthSafetyFactor;
      }

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
		      &JeansLengthRefinement, temperature) == FAIL) {
	    ENZO_FAIL("Error in star_maker3");
      }

    if (STARMAKE_METHOD(INSTANT_STAR)) {

      //---- MODIFIED SF ALGORITHM (NO-JEANS MASS, NO dt DEPENDENCE)

      FORTRAN_NAME(star_maker7)(
       GridDimension, GridDimension+1, GridDimension+2,
       BaryonField[DensNum], dmfield, temperature, BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], cooling_time,
       &dtFixed, BaryonField[NumberOfBaryonFields], BaryonField[MetalNum],
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
          tg->ParticleAttribute[2]);
    } 

      /* Delete any merged particles (Mass == FLOAT_UNDEFINED) */
      
      for (int n = 0; n < NumberOfParticles; n++)
	if (ParticleType[n] == SinkParticleType && 
	    ParticleMass[n] == FLOAT_UNDEFINED)
	  NumberOfStarParticles--;

      this->CleanUpMovedParticles();

    } // ENDIF sinks
 
    /* If not set in the above routine, then set the metal fraction to zero. */
 
    if (MetallicityField == FALSE || StarParticleCreation < 2)
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
	printf("StarParticle: New StarParticles = %"ISYM"\n", NumberOfNewParticles);
 
      /* Set the particle numbers. */
 
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

    int oldDebug = debug;
    if (debug) debug=0;
    delete tg; // temporary grid
    debug = oldDebug;

    //    if (debug) printf("StarParticle: end\n");
 
  }
 

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
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
 
  } // end: if NORMAL_STAR
 
  if (STARFEED_METHOD(UNIGRID_STAR)) {

    //---- UNIGRID (NON-JEANS MASS) VERSION
 
      FORTRAN_NAME(star_feedback3)(
       GridDimension, GridDimension+1, GridDimension+2,
          BaryonField[DensNum], dmfield,
          BaryonField[TENum], BaryonField[GENum], BaryonField[Vel1Num],
          BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[MetalNum],
          BaryonField[MetalNum+1], BaryonField[MetalNum+2],
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
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
 
  } // end: if UNIGRID_STAR
 
  if (STARFEED_METHOD(KRAVTSOV_STAR)) {  

    //---- KRAVTSOV STAR FORMATION ALGORITHM

      FORTRAN_NAME(star_feedback4)(
       GridDimension, GridDimension+1, GridDimension+2,
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
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);

  } // end: if KRAVSTOV STAR

  if (STARFEED_METHOD(INSTANT_STAR)) {

    //---- MODIFIED SF ALGORITHM (NO-JEANS MASS, NO dt DEPENDENCE)
 
      FORTRAN_NAME(star_feedback7)(
       GridDimension, GridDimension+1, GridDimension+2,
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
          ParticleAttribute[2], &RadiationData.IntegratedStarFormation);
 
  } 

  /* Convert the species back from fractional densities to real densities. */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (FieldType[field] >= ElectronDensity && 
	FieldType[field] <= ExtraType1 ) {
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
 
  if (h2field != NULL) delete [] h2field;
  delete [] temperature;
  delete [] dmfield;
  delete [] BaryonField[NumberOfBaryonFields];   // refinement flag field
  BaryonField[NumberOfBaryonFields] = NULL;
 
  //if (debug) printf("StarParticle: end\n");

  JBPERF_STOP("grid_StarParticleHandler");
  return SUCCESS;
}
