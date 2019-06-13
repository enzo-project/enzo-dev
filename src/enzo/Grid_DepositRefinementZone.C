/*********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED IN THE REFINEMENT 
/             ZONE OF AN ACTIVE PARTICLE)
/
/  written by: Nathan Goldbaum
/  date:       January, 2012
/  modified1: Stephen Skory, Sept 2012, name change.
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

FLOAT calc_dist2(FLOAT x1, FLOAT y1, FLOAT z1,
    FLOAT x2, FLOAT y2, FLOAT z2, FLOAT width[]);

int grid::DepositRefinementZone(int level, FLOAT* ParticlePosition,
    FLOAT RefinementRadius)
{
  /* Return if this grid is not on this processor. */

  int dim, method = 0, ParticleMassMethod, i, j, k, NumberOfFlaggedCells = 0, size=1;
  int a;
  float MustRefineMass;
  bool overlaps[GridRank];
  FLOAT CellSize, LeftCorner[MAX_DIMENSION], RightCorner[MAX_DIMENSION];
  FLOAT period[MAX_DIMENSION], left, right, pleft, pright, dist2, rad2;
  
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check whether accretion zone overlaps with the grid Need to
     include ghost zones because the ParticleMassFlaggingField covers
     the ghost zones as well. We need to check for periodicity.*/

  CellSize = CellWidth[0][0];

  rad2 = (RefinementRadius + CellSize) * (RefinementRadius + CellSize);

  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    LeftCorner[dim] = CellLeftEdge[dim][0];
    RightCorner[dim] = LeftCorner[dim] + CellSize*FLOAT((GridDimension[dim]+1));
    period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    overlaps[dim] = true;
  }

  for (dim = 0; dim < GridRank; dim++) {
    left = ParticlePosition[dim] - RefinementRadius;
    right = ParticlePosition[dim] + RefinementRadius;
    // we use or (||) below because and (&&) might exclude smallish grids.
    if ((LeftCorner[dim] > left) || (RightCorner[dim] < right)) {
      overlaps[dim] = false;
    }
    if (overlaps[dim] == false) {
      // check for periodicity.
      pleft = max(period[dim] - fabs(left), left);
      pright = min(fabs(right - period[dim]), right);
      if ((pleft != left) && (RightCorner[dim] > left)) {
        overlaps[dim] = true;
      }
      if ((pright != right) && (LeftCorner[dim] < right)) {
        overlaps[dim] = true;
      }
    }
  } // dim

  for (dim = 0; dim < GridRank; dim++) {
    if (overlaps[dim] == false) return SUCCESS;
  }

  /* Error checks */

  if (ParticleMassFlaggingField == NULL)
    ENZO_FAIL("Particle Mass Flagging Field is undefined!");

  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    if (CellFlaggingMethod[method] == 4)
      ParticleMassMethod = method;
  }

  /* Find mass that will trigger refinement */
  
  MustRefineMass = 1.001*MinimumMassForRefinement[ParticleMassMethod] *
    POW(RefineBy, level * MinimumMassForRefinementLevelExponent[ParticleMassMethod]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the minimimum needed to trigger refinemtn */

  FlaggingField = new int[size];
  for (i = 0; i < size; i++)
    FlaggingField[i] = 0;

  /* Loop over all cells and flag the ones that overlap the accretion zone */

  int index;

  for (i = 0; i < GridDimension[0]; i++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (k = 0; k < GridDimension[2]; k++) {
	index = (k*GridDimension[1]+j)*GridDimension[0]+i;
	dist2 = calc_dist2(CellLeftEdge[0][i] + CellSize/2.,
		  CellLeftEdge[1][j] + CellSize/2.,
		  CellLeftEdge[2][k] + CellSize/2.,
		  ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
		  period);
	if (dist2 <= rad2) {
	  FlaggingField[index] = 1;
	  NumberOfFlaggedCells++;
	}
      } // k
    } // j
  } // i

  /* Set ParticleMassFlaggingField appropriately */

  for (i = 0; i < size; i++)
    ParticleMassFlaggingField[i] += (FlaggingField[i] > 0) ? MustRefineMass : 0;

  delete [] FlaggingField;
  FlaggingField = NULL;

  return SUCCESS;
}
