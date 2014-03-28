/***********************************************************************
/
/  GRID CLASS (DEPOSIT MUST-REFINE PARTICLES)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  May, 2009 by John Wise
/                Works with local particle storage in RebuildHierarchy,
/                which appropriately distributes memory usage.
/  modified2:  Sep, 2009 by Ji-hoon Kim
/                Now PARTICLE_TYPE_MBH also has the MUST_REFINE feature
/
/  PURPOSE:
/
************************************************************************/
 
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
#include "CosmologyParameters.h"

 
/* function prototypes */
 
extern "C" void PFORTRAN_NAME(cic_flag)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, float *partmass, int *ndim, int *npositions,
                        int *itype, int *ffield, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize,
			int *imatch1, int *imatch2, float *minmassmust, int *buffersize,
			float *unipartmass);
 
 
int grid::DepositMustRefineParticles(int pmethod, int level)
{
  /* declarations */
  //printf("grid::DepositMustRefineParticles called \n");
  int i, dim, size = 1, ParticleTypeToMatch1, ParticleTypeToMatch2;
  FLOAT LeftEdge[MAX_DIMENSION], CellSize;
  int ParticleBufferSize;

  ParticleBufferSize = 1;
  if (ProblemType == 106 || ProblemType ==107)
    ParticleBufferSize = 16;


  /* error check */
 
  if (ParticleMassFlaggingField == NULL) {
    fprintf(stderr, "Particle Mass Flagging Field is undefined.\n");
    return -1;
  }

  /* If refining region before supernova (to be safe in its last 5% of
     the lifetime), temporarily set particle type of star to
     PARTICLE_TYPE_MUST_REFINE. */

  if (PopIIISupernovaMustRefine == TRUE)
    this->ChangeParticleTypeBeforeSN(PARTICLE_TYPE_MUST_REFINE, level,
				     &ParticleBufferSize);

  /* Set Left edge of grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }

  CellSize = float(CellWidth[0][0]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the refinement criteron. */

  FlaggingField = new int[size];
  for (i = 0; i < size; i++)
    FlaggingField[i] = 0;

  /* Loop over all the particles, using only particles marked as
     must-refine particles. */

  ParticleTypeToMatch1 = PARTICLE_TYPE_MUST_REFINE;
  ParticleTypeToMatch2 = PARTICLE_TYPE_MBH;
 
  float UniformParticleMass = 0.0;
  if (ProblemType == 30 and MustRefineParticlesCreateParticles == 3){
    float OmegaCDMNow = 1.0 - OmegaLambdaNow;
    UniformParticleMass = OmegaCDMNow/OmegaMatterNow;
  }
  PFORTRAN_NAME(cic_flag)(
	   ParticlePosition[0], ParticlePosition[1], ParticlePosition[2], ParticleMass,
	   &GridRank, &NumberOfParticles, ParticleType, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleTypeToMatch1, &ParticleTypeToMatch2, 
	   &MustRefineParticlesMinimumMass, &ParticleBufferSize,
			  &UniformParticleMass);

  /* Increase particle mass flagging field for definite refinement */

  float MustRefineMass = 
    1.001*MinimumMassForRefinement[pmethod] * 
    POW(RefineBy, level * MinimumMassForRefinementLevelExponent[pmethod]);
  if (ProblemType == 28)
    MustRefineMass = 0;
  
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i]) {
      ParticleMassFlaggingField[i] = MustRefineMass;
      NumberOfFlaggedCells++;
    }

  /* If refining region before supernova, change particle type back to
     its original value. */

  if (PopIIISupernovaMustRefine == TRUE)
    this->ChangeParticleTypeBeforeSN(PARTICLE_TYPE_RESET, level);

  /* Clean up */

  delete [] FlaggingField;
  FlaggingField = NULL;

  return NumberOfFlaggedCells;
 
}
