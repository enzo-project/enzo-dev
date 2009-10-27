/***********************************************************************
/
/  GRID CLASS (SPLIT PARTICLES INTO CHILDREN PARTICLES)
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:  

/  PURPOSE: This routine splits particles into 13 (=12+1) children particles 
/           when requested.  See Kitsionas & Whitworth (2002) for the
/           technical details which was already implemented and used 
/           in the partilce splitting technique of the SPH scheme.
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

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
extern "C" void FORTRAN_NAME(particle_splitter)(int *nx, int *ny, int *nz,
             int *idual, int *imetal, hydro_method *imethod, float *dt, 
		       float *dx, FLOAT *t, float *z, 
             float *d1, float *x1, float *v1, float *t1,
	     FLOAT *xstart, FLOAT *ystart, FLOAT *zstart, int *ibuff,
             int *npart,
             FLOAT *xpold, FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, float *wpold,
	     float *mpold, float *tdpold, float *tcpold, float *metalfold, int *typeold,
	     int *nmax, int *npartnew, int *children, int *level,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
	     float *mp, float *tdp, float *tcp, float *metalf, int *type, int *ran1_init); 

  
int grid::ParticleSplitter(int level)
{

  if (ParticleSplitterIterations == 0)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  /* initialize */
 
  int dim, i, j, k, index, size, field, GhostZones = DEFAULT_GHOST_ZONES;
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
 
  if (MultiSpecies > 1) {
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    H2IINum  = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
  }

  /* Find metallicity field and set flag. */
 
  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum; 
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    ENZO_FAIL("");
  }

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
  int MaximumNumberOfNewParticles = (ChildrenPerParent+1)*NumberOfParticles;
  tg->AllocateNewParticles(MaximumNumberOfNewParticles);

  int NumberOfNewParticles = 0;

  /* ------------------------------------------------------------------- */
  /*                       NOW SPLIT THE PARTICLES                       */
  /* ------------------------------------------------------------------- */

  if(debug)
    fprintf(stdout, "grid::ParticleSplitter:  NumberOfParticles before splitting = %d\n", 
	    NumberOfParticles);
 
  if (NumberOfParticles > 0)
    FORTRAN_NAME(particle_splitter)(
       GridDimension, GridDimension+1, GridDimension+2,
       &DualEnergyFormalism, &MetallicityField, &HydroMethod,
       &dtFixed, &CellWidthTemp,
          &Time, &zred,
       &DensityUnits, &LengthUnits, &VelocityUnits, &TimeUnits,
          CellLeftEdge[0], CellLeftEdge[1], CellLeftEdge[2], &GhostZones,
       &NumberOfParticles, 
       ParticlePosition[0], ParticlePosition[1],
          ParticlePosition[2],
       ParticleVelocity[0], ParticleVelocity[1],
          ParticleVelocity[2],
       ParticleMass, ParticleAttribute[1], ParticleAttribute[0],
       ParticleAttribute[2], ParticleType,
       &MaximumNumberOfNewParticles, &NumberOfNewParticles, &ChildrenPerParent, &level,
       tg->ParticlePosition[0], tg->ParticlePosition[1],
          tg->ParticlePosition[2],
       tg->ParticleVelocity[0], tg->ParticleVelocity[1],
          tg->ParticleVelocity[2],
       tg->ParticleMass, tg->ParticleAttribute[1], tg->ParticleAttribute[0],
       tg->ParticleAttribute[2], tg->ParticleType, &ran1_init);

  /* If not set in the above routine, then set the metal fraction to zero. */
  
  if (MetallicityField == FALSE)
    for (i = 0; i < NumberOfNewParticles; i++)
      tg->ParticleAttribute[2][i] = 0.0;

  /* Set the particle numbers (=indices). Following the convention in Grid_StarParticleHandler, 
     particles won't get *meaningful* indices here;  instead this will be done in 
     CommunicationUpdateStarParticleCount in ParticleSplitter. */
    
  for (i = 0; i < NumberOfNewParticles; i++)
    tg->ParticleNumber[i] = INT_UNDEFINED;

  /* Move any new particles into their new homes. */
 
  if (NumberOfNewParticles > 0) {
    
    if (debug) 
      fprintf(stdout, "grid::ParticleSplitter:  NumberOfNewParticles = %"ISYM"\n", NumberOfNewParticles);    

    /* Move Particles into this grid (set cell size) using the fake grid. */
    
    tg->NumberOfParticles = NumberOfNewParticles;
    for (dim = 0; dim < GridRank; dim++) {
      tg->CellWidth[dim] = new FLOAT[1];
      tg->CellWidth[dim][0] = CellWidth[dim][0];
    }
    this->MoveAllParticles(1, &tg);

    if (debug)
      fprintf(stdout, "grid::ParticleSplitter:  NumberOfParticles (New) = %"ISYM"\n", NumberOfParticles);    
    
  } // end: if (NumberOfNewParticles > 0)


  /* Clean up. */
  
  delete tg; // temporary grid
  delete [] BaryonField[NumberOfBaryonFields];   // refinement flag field
  BaryonField[NumberOfBaryonFields] = NULL;

  LCAPERF_STOP("grid_ParticleSplitter");
  return SUCCESS;

}
