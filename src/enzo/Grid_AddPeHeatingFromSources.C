/***********************************************************************
/
/  ADD PE HEATING RATE FROM FUV RADIATION FROM SHINING PARTICLES
/
/  written by: Andrew Emerick
/  date:       December, 2016
/       taken from Grid_AddH2DissociationFromSources
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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int f, int farray[], int n);

float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e, const float &Z,
                                      const float &T, const float &G);

int grid::AddPeHeatingFromSources(Star *AllStars)
{

  Star *cstar;
  FLOAT DomainWidth[MAX_DIMENSION];
  FLOAT *ddr2[MAX_DIMENSION];
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS];
  int ipart, dim, i, j, k, index, indixe, nbins;
  int ActiveDims[MAX_DIMENSION];
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  const double pc = 3.086e18, clight = 3e10;
  const double eV_erg = 6.241509e11;
  const double m_e = 9.109E-28; // in g
  const double m_h = 1.673E-24; // in g


  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddPeHeating");

  /* Find Multi-species fields. */

  if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
                                  HeIIINum, HMNum, H2INum, H2IINum, DINum, 
                                  DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  /* Get photo-ionization fields */

  int PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);

  /* Initialize the Pe heating rate field */
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (AllStars == NULL)
    return SUCCESS;

  if(ProblemType != 50 && !(STARMAKE_METHOD(INDIVIDUAL_STAR)))
    return SUCCESS;


  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1, EnergyUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, PhotonTime);
  EnergyUnits = DensityUnits * VelocityUnits * VelocityUnits;

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  /* get temperature field */
  float *temperature;

  temperature = new float[size];

  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
  }

  for (dim = 0; dim < GridRank; dim++){
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ddr2[dim] = new FLOAT[ActiveDims[dim]];
  }

  // Dilution factor to prevent breaking of rate solver near the star
  float dilutionRadius = 0.125 * this->CellWidth[0][0];
  float dilRadius2     = dilutionRadius * dilutionRadius;
  float LightTravelDist = TimeUnits * clight / LengthUnits;

  float PeConversion = 1.0 / (EnergyUnits / TimeUnits);
  float FUVLuminosity = 0.0;

  if (ProblemType == 50) ENZO_FAIL("Ptype = 50 not implemented in PeHeating");

  for (cstar = AllStars; cstar; cstar = cstar->NextStar){

    // Skip if not 'living'
    // checks are different for individual star particles
    if (STARMAKE_METHOD(INDIVIDUAL_STAR)){
      // these checks shouldn't be needed since things are
      // zeroed, but saves some time rather than running
      // through loops below

      if ( (cstar->type != PARTICLE_TYPE_INDIVIDUAL_STAR) ||
           (cstar->BirthMass < IndividualStarOTRadiationMass ))
       continue;

    } else {
      if (!(cstar->FeedbackFlag == NO_FEEDBACK ||
          cstar->FeedbackFlag == CONT_SUPERNOVA))
      continue;
    }

    /* Determine FUV rates */
    if (cstar->ComputePhotonRates(TimeUnits, nbins, energies, Luminosity) == FAIL){
      ENZO_FAIL("Error in ComputePhotonRates from AddPeHeatingFromSources.\n");
    }
    /* this->Luminosity is photon / s, energies is in eV */
    FUVLuminosity = (Luminosity[4]*energies[4]) / (4.0 * M_PI * eV_erg);

    /* Pre-calculate distances from cells to source */
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
           i++, index++) {

        // Calculate dr_i first, then square it
        ddr2[dim][i] = 
          fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][index] -
               cstar->pos[dim]);
//        ddr2[dim][i] = m1in(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
        ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
      }

   /* Loop over cells */

    double radius2, radius2_yz;
    double FUVflux;
    for (k = 0; k < ActiveDims[2]; k++) {
      for (j = 0; j < ActiveDims[1]; j++) {
        radius2_yz = ddr2[1][j] + ddr2[2][k];
        index = GRIDINDEX(0, j, k);
        for (i = 0; i < ActiveDims[0]; i++, index++) {
          radius2 = radius2_yz + ddr2[0][i];

          if (radius2 < dilRadius2) // need r^2 in cgs
            FUVflux = FUVLuminosity / (dilRadius2 * LengthUnits * LengthUnits);
          else
            FUVflux = FUVLuminosity / (radius2 * LengthUnits * LengthUnits);

          float n_H, n_e, Z;

          n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

          if ( MultiSpecies > 1){ /* include H2 */
             n_H += this->BaryonField[HMNum][index] +
                      0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
          }


          n_H *= DensityUnits / m_h;

          n_e  = this->BaryonField[DeNum][index] * DensityUnits / m_e;

          Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens


          BaryonField[PeNum][index] += ComputeHeatingRateFromDustModel(n_H, n_e, Z,
                                                                       temperature[index], FUVflux) * PeConversion;
            //} // ENDIF
        } // END: i-direction
      } // END: j-direction
    } // END: k-direction
  } // ENDFOR stars

  for (dim = 0; dim < GridRank; dim++)
    delete [] ddr2[dim];

  delete [] temperature;

  return SUCCESS;

}
