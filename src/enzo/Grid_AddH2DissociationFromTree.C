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
#include "phys_constants.h"

#define MIN_OPENING_ANGLE 0.2  // 0.2 = arctan(11.3 deg)

// Defined in Grid_AddH2DissociationFromSources
static double CalculateH2IICrossSection(float Energy);
static double CalculateIRCrossSection(float Energy);

// Defined in FindSuperSourceByPosition.C:
float CalculateLWFromTree(const FLOAT pos[], const float angle,
        const SuperSourceEntry *Leaf, const float min_radius,
        float result0);

float CalculateIRFromTree(const FLOAT pos[], const float angle,
        const SuperSourceEntry *Leaf, const float min_radius,
        float result0);

float CalculateFUVFromTree(const FLOAT pos[], const float angle,
        const SuperSourceEntry *Leaf, const float min_radius,
        float result0);

int FindSuperSourceByPosition(FLOAT *pos, SuperSourceEntry **result,
			      int DEBUG);
//---------------------------------------------------------

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddH2DissociationFromTree(void)
{

  int i, j, k, index, dim, ci;
  FLOAT pos[MAX_DIMENSION];
  FLOAT radius2;
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  float kdiss_r2;
  double LWLuminosity, H2ISigma = 3.71e-18;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if there are no sources */

  if (SourceClusteringTree == NULL)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddH2Dissociation");

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum,
				  kdissH2INum, kphHMNum, kdissH2IINum);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits,
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 10.0 * pc_cm / (double) LengthUnits;
  float dilRadius2 = dilutionRadius * dilutionRadius;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);
  double LConv_inv = 1.0 / LConv;

  /* Find sources in the tree that contribute to the cells */

  SuperSourceEntry *Leaf;
  const float factor    = LConv_inv / (4.0 * pi);
  const float H2Ifactor = factor * H2ISigma;
  float angle;

  Leaf = SourceClusteringTree;


  // compute cross sections
  //   assuming all energies are the same for each band always
  //   this is a strong assumption and this should be fixed to take the actual
  //   energies from RadiationSourceEntry
  float H2IICrossSection[3];
  float HMCrossSection[3];
  float PhotonEnergy[3] = {IR_photon_energy, FUV_photon_energy, LW_photon_energy};

  // all cross sections include the above RT conversion factor for simplicity
  // in looping over sources below
  for (int n = 0; n < 3; n++){
    H2IICrossSection[n] = factor * CalculateH2IICrossSection(PhotonEnergy[n]);
    HMCrossSection[n]   = factor * CalculateIRCrossSection(PhotonEnergy[n]);
  }

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

        // Compute IR, FUV, and LW bands as all three can affect H2 or HM
        double IRLuminosity  = CalculateIRFromTree (pos, angle, Leaf, dilRadius2, 0);
        double FUVLuminosity = CalculateFUVFromTree(pos, angle, Leaf, dilRadius2, 0);
        double LWLuminosity  = CalculateLWFromTree (pos, angle, Leaf, dilRadius2, 0);

        // H2I dissociation only from LW radiationn
        //   this one is simple
        BaryonField[kdissH2INum][index] = LWLuminosity * H2Ifactor;

        // H2II and HM dissociation can occur from multiple bands.
        // need to compute this for each band:
        //     AJE:    hard coding energies here is generally bad
        //             this really should be set up to be gaurunteed to
        //             match RadiationSourceEntry


        // H2II dissociation from LW, FUV, and IR
        BaryonField[kdissH2IINum][index] = IRLuminosity  * H2IICrossSection[0] +
                                           FUVLuminosity * H2IICrossSection[1] +
                                           LWLuminosity  * H2IICrossSection[2];

        // HM dissociation from FUV and IR
        BaryonField[kphHMNum][index]     = IRLuminosity  * HMCrossSection[0] +
                                           FUVLuminosity * HMCrossSection[1];

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return SUCCESS;

}


static double CalculateH2IICrossSection(float Energy)
{
  float X = 0.0;
  float a = 3.35485518, b = 0.93891875, c = -0.01176537;
  /* Fit created from Stancil et al. 1994 */
  X = log(Energy/8.0);
  return a*pow(10.0, -b*X*X)*exp(-c*X) * 1e-18;

}
static double CalculateIRCrossSection(float Energy)
{
  float X = 0.0;
  float A = 3.486e-16;
  /* Fit taken from Tegmark et al. (1997) */
  X = Energy / 0.74;
  return (A*pow((X - 1), 1.5)/pow(X, 3.11));
}
