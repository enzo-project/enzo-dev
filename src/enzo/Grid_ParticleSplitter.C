/***********************************************************************
/
/  GRID CLASS (SPLIT PARTICLES INTO CHILDREN PARTICLES)
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:  

/  PURPOSE: This routine splits particles into 13 (=12+1) children particles 
/           when requested.  See Kitsionas & Whitworth (2002) for the
/           technical details of particle splitting,  which was already 
/           implemented and used in SPH/Gadget.
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

#define NO_DEBUG_PS 

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
extern "C" void FORTRAN_NAME(particle_splitter)(int *nx, int *ny, int *nz,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
	     float *r, float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
	     FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff,
             int *npart,
             FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, float *wpold,
	     float *mpold, float *tdpold, float *tcpold, float *metalfold, int *typeold,
	     int *nmax, int *npartnew, int *children, int *level,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int *type, 
	     int *iterations, float *separation, int *ran1_init, 
	     FLOAT *rr_leftedge, FLOAT *rr_rightedge); 

  
int grid::ParticleSplitter(int level)
{

  if (ParticleSplitterIterations == 0)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  if (GridRank <=2)
    ENZO_FAIL("GridRank <= 2 has never been tested; do you really want to continue?");

  /* Initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = NumberOfGhostZones;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,H2INum, H2IINum;

  LCAPERF_START("grid_ParticleSplitter");
 
  /* Compute size (in floats) of the current grid. */
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  this->DebugCheck("ParticleSplitter");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
 
  /* Find metallicity field and set flag. */
 
  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
              MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;

  /* Compute the redshift. */
 
  float zred;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  zred = 1.0*(1.0+InitialRedshift)/a - 1.0;
 
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
 
  float CellWidthTemp = float(CellWidth[0][0]);

  /* Generate a fake grid to keep the particles in. */
  
  grid *tg = new grid;
  tg->GridRank = GridRank;
  tg->ProcessorNumber = ProcessorNumber;
  
  /* Allocate space for new particles. */
  
  int ChildrenPerParent = 12;  //12+1 = 13 will be the final number of particles per parent
  int MaximumNumberOfNewParticles = ChildrenPerParent*NumberOfParticles+1;
  tg->AllocateNewParticles(MaximumNumberOfNewParticles);

  int NumberOfNewParticles = 0;

  /* ------------------------------------------------------------------- */
  /*                       NOW SPLIT THE PARTICLES                       */
  /* ------------------------------------------------------------------- */

#ifdef DEBUG_PS
  fprintf(stdout, "grid::ParticleSplitter:  NumberOfParticles before splitting = %d, MyProcessorNumber = %d\n", 
	  NumberOfParticles, MyProcessorNumber); 
#endif
 
  if (NumberOfParticles > 0) {

#define NO_PARTICLE_IN_GRID_CHECK 

#ifdef PARTICLE_IN_GRID_CHECK
    int xindex, yindex, zindex;
    for (i = 0; i < NumberOfParticles; i++) {
      
      xindex = (int)((ParticlePosition[0][i] - CellLeftEdge[0][0]) / CellWidthTemp);
      yindex = (int)((ParticlePosition[1][i] - CellLeftEdge[1][0]) / CellWidthTemp); 
      zindex = (int)((ParticlePosition[2][i] - CellLeftEdge[2][0]) / CellWidthTemp); 

      if (xindex < 0 || xindex >= GridDimension[0] || 
	  yindex < 0 || yindex >= GridDimension[1] || 
	  zindex < 0 || zindex >= GridDimension[2])
	fprintf(stdout, "grid::PS: parent particle out of grid (C level): xind, yind, zind, level = %d, %d, %d, %d\n",
		xindex, yindex, zindex, level); 
    }
#endif
    if(!tg->CreateChildParticles(CellWidthTemp, NumberOfParticles, ParticleMass,
				 ParticleType, ParticlePosition, ParticleVelocity,
				 ParticleAttribute, CellLeftEdge, GridDimension, 
				 MaximumNumberOfNewParticles, &NumberOfNewParticles))
      {
	fprintf(stdout, "Failed to create child particles in grid %d\n", this->GetGridID());
	return FAIL;
      }

  }

  /* Move any new particles into their new homes. */
 
  if (NumberOfNewParticles > 0) {

#ifdef DEBUG_PS
    fprintf(stdout, "grid::ParticleSplitter:  NumberOfNewParticles = %d, MyProcessorNumber = %d\n", 
	    NumberOfNewParticles, MyProcessorNumber);    
#endif

    /* If not set in the above routine, then set the metal fraction to zero. */
    
    if (MetallicityField == FALSE && NumberOfParticleAttributes > 1)
      for (i = 0; i < NumberOfNewParticles; i++)
	tg->ParticleAttribute[2][i] = 0.0;
    
    /* Set the particle numbers (=indices). Following the convention 
       in Grid_StarParticleHandler, particles won't get indices here;  
       instead it will be done in CommunicationUpdateStarParticleCount in ParticleSplitter. 
       Plus, because we create different types of particles, 
       here we set ParticleNumber differently so that they can be dealt with later. */
    
    for (i = 0; i < NumberOfNewParticles; i++)
	tg->ParticleNumber[i] = INT_UNDEFINED;
    
    /* Move Particles into this grid (set cell size) using the fake grid. */
    
    tg->NumberOfParticles = NumberOfNewParticles;
    for (dim = 0; dim < GridRank; dim++) {
      tg->CellWidth[dim] = new FLOAT[1];
      tg->CellWidth[dim][0] = CellWidth[dim][0];
    }
    this->MoveAllParticles(1, &tg);

#ifdef DEBUG_PS
    fprintf(stdout, "grid::ParticleSplitter:  NumberOfParticles(New) = %d, MyProcessorNumber = %d\n", 
	    NumberOfParticles, MyProcessorNumber);    
#endif
    
  } // end: if (NumberOfNewParticles > 0)


  /* Clean up. */
  
  delete tg; // temporary grid

  LCAPERF_STOP("grid_ParticleSplitter");
  return SUCCESS;

}
