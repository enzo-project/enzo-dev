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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int f, int farray[], int n);

float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      // const float &T,
                                      const float &Z, const float &G, const float &dx);

int grid::AddPeHeatingFromSources(Star *AllStars)
{

  Star *cstar;
  FLOAT DomainWidth[MAX_DIMENSION];
  FLOAT *ddr2[MAX_DIMENSION];
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS];
  int dim, i, j, k, index, nbins;
  int ActiveDims[MAX_DIMENSION];
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;


  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddPeHeating");

  /* Find Multi-species fields. */
// done in RT initialize  this->ZeroPhotoelectricHeatingField();

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
  int FUVRateNum = -1;
  FUVRateNum = FindField(FUVRate, this->FieldType, this->NumberOfBaryonFields);


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
  const float dilutionRadius = 0.25 * this->CellWidth[0][0];
  const float dilRadius2     = dilutionRadius * dilutionRadius;

  const float PeConversion = 1.0 / (EnergyUnits / TimeUnits);
  float FUVLuminosity = 0.0;
  // Note inconsistency here (same as defined elsewhere, but different form.
  const float FluxConv = EnergyUnits / TimeUnits * LengthUnits;
  const float FluxConv_inv = 1.0 / FluxConv;


  if (ProblemType == 50) ENZO_FAIL("Ptype = 50 not implemented in PeHeating");

  const double clight_code = clight * TimeUnits / LengthUnits;

  for (cstar = AllStars; cstar; cstar = cstar->NextStar){

    // Skip if not 'living'
    // checks are different for individual star particles
    if (STARMAKE_METHOD(INDIVIDUAL_STAR)){
      // these checks shouldn't be needed since things are
      // zeroed, but saves some time rather than running
      // through loops below

      if ( !((cstar->type == PARTICLE_TYPE_INDIVIDUAL_STAR) ||
             (cstar->type == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII)) ||
             (cstar->BirthMass < IndividualStarOTRadiationMass) )
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
    // compute FUV uminosity gives FUV luminosity in erg / s

    if (abs(cstar->type) != IndividualStarPopIII){
      cstar->ComputeFUVLuminosity(FUVLuminosity);
    } else {
      FUVLuminosity = Luminosity[7]*energies[7] / eV_erg; // now in erg/s
    }

    /* Pre-calculate distances from cells to source */
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim];
           i++, index++) {

        // Calculate dr_i first, then square it
        ddr2[dim][i] =
          fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][index] -
               cstar->pos[dim]);
        if (RadiativeTransferPeriodicBoundary) ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
        ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
      }

   /* Loop over cells */

    double radius2, radius2_yz;
    double FUVflux = FUVLuminosity / (4.0 * pi * LengthUnits * LengthUnits); // mostly converted to flux
    for (k = 0; k < ActiveDims[2]; k++) {
      for (j = 0; j < ActiveDims[1]; j++) {
        radius2_yz = ddr2[1][j] + ddr2[2][k];
        index = GRIDINDEX(0, j, k);
        for (i = 0; i < ActiveDims[0]; i++, index++) {
          radius2 = radius2_yz + ddr2[0][i];

          float max_distance = (this->Time - cstar->ReturnBirthTime()) * clight_code;

          if ( sqrt(radius2) > max_distance) continue; // does not contribute

          double LocalFUVflux = 0.0;
          if (radius2 < dilRadius2){ // need r^2 in cgs - this is done above in FUVflux
            LocalFUVflux = FUVflux / (dilRadius2);
          } else{
            LocalFUVflux = FUVflux / (radius2);
          }
          // LocalFUVFlux now in cgs units (erg/s/cm^2)

          float n_H, n_e, Z;

          n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

          if ( MultiSpecies > 1){ /* include H2 */
             n_H += this->BaryonField[HMNum][index] +
                      0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
          }


          n_H *= DensityUnits / mh;

          n_e  = this->BaryonField[DeNum][index] * DensityUnits / me;

          Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens

          if (FUVRateNum > 0){
            BaryonField[FUVRateNum][index] += LocalFUVflux * FluxConv_inv;
          }


          if (temperature[index] > IndividualStarFUVTemperatureCutoff) {
            BaryonField[PeNum][index] = 0.0;
          } else {
            /* Need to make decision on whether or not to place T cut here or in Grackle wrapper */

            BaryonField[PeNum][index] += ComputeHeatingRateFromDustModel(n_H, n_e,
                                                               // temperature[index],
                                                                       Z, LocalFUVflux,
                                                                       this->CellWidth[0][0]*LengthUnits) * PeConversion;
          } // end T check

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
