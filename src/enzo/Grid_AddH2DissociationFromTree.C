/***********************************************************************
/
/  ADD H2 DISSOCIATION EMISSION FROM STAR PARTICLES FROM A TREE
/
/  written by: John Wise
/  date:       March, 2011
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
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
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "Star.h"

#define MIN_OPENING_ANGLE 0.2  // 0.2 = arctan(11.3 deg)

float CalculateLWFromTree(const FLOAT pos[], const float angle, 
			  const SuperSourceEntry *Leaf, const float min_radius, 
			  float result0);
int FindSuperSourceByPosition(FLOAT *pos, SuperSourceEntry **result,
			      int DEBUG);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddH2DissociationFromTree(void)
{

  const double pc = 3.086e18, clight = 3e10;

  int i, j, k, index, dim, ci;
  FLOAT pos[MAX_DIMENSION];
  FLOAT radius2;
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS], kdiss_r2;
  double H2Luminosity, H2ISigma = 3.71e-18;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if there are no sources */

  if (SourceClusteringTree == NULL)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddH2Dissociation");

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum, 
				  kdissH2INum);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 10.0 * pc / (double) LengthUnits;
  float dilRadius2 = dilutionRadius * dilutionRadius;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);
  double LConv_inv = 1.0 / LConv;

  /* Find sources in the tree that contribute to the cells */

  SuperSourceEntry *Leaf;
  float factor = LConv_inv * H2ISigma / (4.0 * M_PI);
  float angle;

  Leaf = SourceClusteringTree;

  // We want to use the source separation instead of the merging
  // radius.  The leaves store the merging radius (ClusteringRadius),
  // so we multiply the angle by merge radius.
  angle = MIN_OPENING_ANGLE * RadiativeTransferPhotonMergeRadius;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

	/* Find the leaves that have an opening angle smaller than
	   the specified minimum and only include those in the
	   calculation */

	H2Luminosity = CalculateLWFromTree(pos, angle, Leaf, dilRadius2, 0);
	BaryonField[kdissH2INum][index] = H2Luminosity * factor;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}
