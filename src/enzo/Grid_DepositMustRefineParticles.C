/***********************************************************************
/
/  GRID CLASS (DEPOSIT MUST-REFINE PARTICLES)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  May, 2009 by John Wise
/                Works with local particle storage in RebuildHierarchy,
/                which appropriately distributes memory usage.
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
extern "C" void FORTRAN_NAME(cic_flag)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        int *itype, int *ffield, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize,
			int *imatch);
 
 
int grid::DepositMustRefineParticles(int pmethod, int level)
{
  /* declarations */
 
  int i, dim, ParticleTypeToMatch, size = 1;
  FLOAT LeftEdge[MAX_DIMENSION], CellSize;
 
  /* error check */
 
  if (ParticleMassFlaggingField == NULL) {
    fprintf(stderr, "Particle Mass Flagging Field is undefined.\n");
    return -1;
  }

  /* Set Left edge of grid. */
 
  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }
  ParticleTypeToMatch = PARTICLE_TYPE_MUST_REFINE;
  CellSize = float(CellWidth[0][0]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the refinement criteron. */

  FlaggingField = new int[size];
  for (i = 0; i < size; i++)
    FlaggingField[i] = 0;

  /* Loop over all the particles, using only particles marked as
     must-refine particles. */
 
  FORTRAN_NAME(cic_flag)(
           ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, ParticleType, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleTypeToMatch);

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

 
  /* Clean up */

  delete [] FlaggingField;
  FlaggingField = NULL;

  return NumberOfFlaggedCells;
 
}
