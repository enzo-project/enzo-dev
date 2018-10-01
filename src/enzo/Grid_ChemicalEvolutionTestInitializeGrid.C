/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A TEST OF CHEMICAL EVOLUTION MODELS)
/
/  written by: Andrew Emerick
/  date:       Feb, 2016
/
/  PURPOSE: Plain, uniform gas grid. Star is deposited in individual_star_maker
/           this is a boring set up.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"


int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int grid::ChemicalEvolutionTestInitializeGrid(float GasDensity, float GasTemperature,
                                              float GasMetallicity){

  if (ProcessorNumber != MyProcessorNumber){
    return SUCCESS;
  }


  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                              &TimeUnits, &VelocityUnits, Time) == FAIL){
      ENZO_FAIL("Error in GetUnits.");
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL){
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum, PeNum, OTLWkdissH2INum;

  if(MultiSpecies){
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                              HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL){
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
  }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1){
    MetallicityField = TRUE;
  } else {
    MetalNum = 0;
    printf("ChemicalEvolutionTest: Metallicity Field not found.\n");
  }

  if(STARMAKE_METHOD(INDIVIDUAL_STAR) && IndividualStarFUVHeating){
    PeNum = FindField(PeHeatingRate, FieldType, NumberOfBaryonFields);
    if (PeNum <= 0){
      ENZO_FAIL("Error identifying pe heating rate field\n");
    }
  }
  if(STARMAKE_METHOD(INDIVIDUAL_STAR) && IndividualStarLWRadiation){
    OTLWkdissH2INum = FindField(OTLWkdissH2I, FieldType, NumberOfBaryonFields);
    if (OTLWkdissH2INum <= 0){
      ENZO_FAIL("Error identifying kdissH2I field\n");
    }
  }

  int size = 1, i;

  for (int dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
  }

  if (TestProblemData.UseMetallicityField){
    for ( i = 0; i < size; i++){
      BaryonField[MetalNum][i] = GasMetallicity * GasDensity;
    }
  }


  if (MultiSpecies){
    // set set background to primordial and 100% ionized (only HII and HeIII)
    for( i = 0; i < size; i++){
      BaryonField[HIINum][i]   = GasDensity * TestProblemData.HydrogenFractionByMass
                                            * TestProblemData.HII_Fraction ;
      BaryonField[HeIINum][i]  = GasDensity * TestProblemData.HeII_Fraction *
                                          (1.0 - TestProblemData.HydrogenFractionByMass);
      BaryonField[HeIIINum][i] = GasDensity * TestProblemData.HeIII_Fraction
                                            * (1.0 - TestProblemData.HydrogenFractionByMass);
      BaryonField[HeINum][i]   = (1.0 - TestProblemData.HydrogenFractionByMass) * GasDensity -
                                 BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];
      if(MultiSpecies > 1){
        BaryonField[HMNum][i]  = TestProblemData.HM_Fraction  * BaryonField[HIINum][i];
        BaryonField[H2INum][i] = TestProblemData.H2I_Fraction * GasDensity
                                                   * TestProblemData.HydrogenFractionByMass;
        BaryonField[H2IINum][i] = TestProblemData.H2II_Fraction * 2.0 * BaryonField[HIINum][i];
      }

      BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass *GasDensity
                                                - BaryonField[HIINum][i];

      if( MultiSpecies > 1){
        BaryonField[HINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i] +
                                 BaryonField[H2INum][i]);
      }

      // Electron density: sum up ionized species
      BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25 * BaryonField[HeIINum][i] +
                                                      0.50 * BaryonField[HeIIINum][i];
      if (MultiSpecies > 1){
        BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i];
      }

      if (MultiSpecies > 2){
        BaryonField[DINum ][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
        BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
        BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    }// loop over cells
  } // Multispecies



  /* Loop over all requested stellar yields species and assign initial values */
  if (TestProblemData.MultiMetals == 2){
    for( int sp = 0; sp < StellarYieldsNumberOfSpecies; sp++){
      if(StellarYieldsAtomicNumbers[sp] > 2){
        int   field_num;
        float fraction;

        fraction = TestProblemData.ChemicalTracerSpecies_Fractions[sp];

        IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[sp]);

        for(i = 0; i < size; i++){ // assign initial values
          BaryonField[field_num][i] = fraction * GasDensity;
        }
      }
    } // loop over yields

  } // MULTI METALS

  if(STARMAKE_METHOD(INDIVIDUAL_STAR) && IndividualStarFUVHeating){
    for(i = 0; i < size; i ++){
      BaryonField[PeNum][i] = 0.0;
    }
  }
  if(STARMAKE_METHOD(INDIVIDUAL_STAR) && IndividualStarLWRadiation){
    for(i = 0; i < size; i ++){
      BaryonField[OTLWkdissH2INum][i] = 0.0;
    }
  }


  /* now go through and initialize the particles */
  int MaximumNumberOfNewParticles = 10000;
  int NumberOfNewParticles = 0;
  this->AllocateNewParticles(MaximumNumberOfNewParticles);

  this->chemical_evolution_test_star_deposit(&MaximumNumberOfNewParticles,
                                             &NumberOfNewParticles, this->ParticleMass,
                                             this->ParticleType, this->ParticlePosition,
                                             this->ParticleVelocity, this->ParticleAttribute);

  if (NumberOfNewParticles > 0) {
    this->NumberOfParticles = NumberOfNewParticles;
  } else if (ChemicalEvolutionTestNumberOfStars > 0) {
    ENZO_FAIL("Was not able to deposit stars in chemical evolution test\n");
  } // end: if (NumberOfNewParticles > 0)



  return SUCCESS;
}
