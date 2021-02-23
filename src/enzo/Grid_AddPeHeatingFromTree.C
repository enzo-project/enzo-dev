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

float CalculateFUVFromTree(const FLOAT pos[], const float angle,
                          const SuperSourceEntry *Leaf, const float min_radius,
                          float result0);
int FindSuperSourceByPosition(FLOAT *pos, SuperSourceEntry **result,
                              int DEBUG);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      // const float &T,
                                      const float &Z,
                                      const float &G, const float &dx);



int grid::AddPeHeatingFromTree(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if there are no sources */

  if (SourceClusteringTree == NULL)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddPeHeating");

  int PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);
  int ElectronNum = FindField(ElectronDensity, this->FieldType, this->NumberOfBaryonFields);

  int FUVRateNum = -1;
  FUVRateNum = FindField(FUVRate, this->FieldType, this->NumberOfBaryonFields);

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }


  /* identify species fields if they exist for proper computation of Mu */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits,
    TimeUnits, aUnits = 1, EnergyUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, PhotonTime);
  EnergyUnits = DensityUnits * VelocityUnits * VelocityUnits;

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* get temperature field */
  float *temperature;
  temperature = new float[size];
  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
  }

  // Dilution factor (prevent breaking in rate solver near star)
  float dilutionRadius = this->CellWidth[0][0] * 0.25;
  float dilRadius2     = dilutionRadius * dilutionRadius;

  /* Find sources in the tree that contribute to the cells */
  SuperSourceEntry *Leaf;
  double LConv = (double) TimeUnits / pow(LengthUnits,3); // this is silly - unconvert a conversion
  double PeConversion = 1.0 / ((double) EnergyUnits / TimeUnits);
  /* Need to include FUV photon energy here since FUVLuminosity in below from source clustering tree
     returns the PHOTON luminosit (1/s) not energy luminosity (erg/s) */
  float factor = (FUV_photon_energy) / ( LConv * eV_erg * (4.0 * pi) *LengthUnits * LengthUnits);
  float angle;
  FLOAT pos[MAX_DIMENSION];

  Leaf = SourceClusteringTree;

  // We want to use the source seperation instead of the merging
  // radius. The leaves store the merging radius (ClusteringRadius)
  // so we multiply the angle by merge radius
  angle = MIN_OPENING_ANGLE * RadiativeTransferPhotonMergeRadius;

  float FUVflux = 0.0;

  // Note inconsistency here (same as defined elsewhere, but different form.
  const float FluxConv = EnergyUnits / TimeUnits * LengthUnits;
  const float FluxConv_inv = 1.0 / FluxConv;


  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++){
    pos[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k];
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++){
      pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++){
        pos[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];

        /* Find the leaves that have an opening angle smaller than
           the specified minimum and only include those in the
           calculation */

        /* FUV from tree should be returning energy flux in RT units,
           so after the conversion factor, this should be Flux in erg / s / cm^2,
           which is what is needed for computing the PE heating rate */
        FUVflux = CalculateFUVFromTree(pos, angle, Leaf, dilRadius2, 0) * factor;

        float n_H, n_e, Z;

        n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

        if ( MultiSpecies > 1){ /* include H2 */
           n_H += this->BaryonField[HMNum][index] +
                    0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
        }


        n_H *= DensityUnits / mh;

        n_e  = this->BaryonField[ElectronNum][index] * DensityUnits / me;

        Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens

        if (FUVRateNum > 0){
          BaryonField[FUVRateNum][index] += (FUVflux * FluxConv_inv);
        }

         /* Need to make decision on whether or not to place T cut here or in Grackle wrapper */

        // assign heating rate from model
        if (temperature[index] > IndividualStarFUVTemperatureCutoff){
          BaryonField[PeNum][index] = 0.0;
        } else {
          BaryonField[PeNum][index]  += ComputeHeatingRateFromDustModel(n_H, n_e,
                                                                       // temperature[index],
                                                                       Z, FUVflux,
                                                                       this->CellWidth[0][0]*LengthUnits) * PeConversion;
        }

      }
    }
  } // end loop

  delete [] temperature;

  return SUCCESS;
}
