/* ---------------------------------------------------------------------
 *
 * Routine to update grid averaged abundances for use with stellar wind
 * mass loading scheme in order to maintain consistent abundances
 *
 * A. Emerick - Sep 2016
 * ---------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"

#include "IndividualStarProperties.h"
#include "StellarYieldsRoutines.h"

#ifndef MU_METAL
# define MU_METAL 16.0
#endif

int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits,
             float *TimeUnits, float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

int UpdateAveragedAbundances(TopGridData *MetaData,
                             LevelHierarchyEntry *LevelArray[],
                             int level, Star* &AllStars){

  if (AllStars == NULL)
    return SUCCESS;

  if ( !IndividualStarStellarWinds ||
       !IndividualStarUseWindMixingModel)
    return SUCCESS;

  LevelHierarchyEntry *Temp;

//  LCAPERF_START("UpdateAveragedAbundances");

  for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l++){

    for( Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel){

      if(Temp->GridData->isLocal()){
        Temp->GridData->CalculateAverageAbundances();
      }
    }
  } // end loop over levels
//  LCAPERF_END("UpdateAveragedAbundances");
  return SUCCESS;
}

int grid::CalculateAverageAbundances(void){
/*

 Loop over the grid, calculating the average abundance
 for each species. This is used for modelling shell mass loading
 in stellar wind ejecta


*/

  float mass_counter[StellarYieldsNumberOfSpecies + 1]; // metal + tracers
  float mass_counter_hot[StellarYieldsNumberOfSpecies + 1];
  float total_mass = 0.0, total_mass_hot = 0.0;
  float *temperature;

  for(int im = 0 ; im < StellarYieldsNumberOfSpecies + 1; im++){
    mass_counter[im] = 0.0;
    mass_counter_hot[im] = 0.0;
  }

  /* get multispecies fields */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, ElectronNum, DensNum, MetalNum;
  if (MultiSpecies){
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  // find fields for density, metal density, electron density, and the heating rate
  DensNum     = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
  MetalNum    = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  int size = this->GridDimension[0] * this->GridDimension[1] * this->GridDimension[2];

  temperature = new float[size];

  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from UpdateAveragedAbundances");
  }

  for(int index = 0; index < size; index++){
    int field_num;

    if( temperature[index] < IndividualStarWindTemperature){
      mass_counter[0] += BaryonField[MetalNum][index];
      for(int im = 0; im < StellarYieldsNumberOfSpecies; im++){

        switch(StellarYieldsAtomicNumbers[im]){

          case 1: // Hydrogen
            mass_counter[im+1] += BaryonField[HINum][index] + BaryonField[HIINum][index];
            if (MultiSpecies > 1){
              mass_counter[im+1] += BaryonField[HMNum][index] + BaryonField[H2INum][index] +
                                    BaryonField[H2IINum][index];
            }
            break;

          case 2: // Helium
            mass_counter[im+1] += BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                                  BaryonField[HeIIINum][index];
            break;

          default:

            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[im]);
            mass_counter[im+1] += BaryonField[field_num][index];

        } // end switch
      } // end species loop
      total_mass += BaryonField[DensNum][index];

    } else {
      mass_counter_hot[0] += BaryonField[MetalNum][index];
      for(int im = 0; im < StellarYieldsNumberOfSpecies; im++){

        switch(StellarYieldsAtomicNumbers[im]){

          case 1: // Hydrogen
            mass_counter_hot[im+1] += BaryonField[HINum][index] + BaryonField[HIINum][index];
            if (MultiSpecies > 1){
              mass_counter_hot[im+1] += BaryonField[HMNum][index] + BaryonField[H2INum][index] +
                                    BaryonField[H2IINum][index];
            }
            break;

          case 2: // Helium
            mass_counter_hot[im+1] += BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                                  BaryonField[HeIIINum][index];
            break;

          default:

            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[im]);
            mass_counter_hot[im+1] += BaryonField[field_num][index];

        } // end switch
      } // end species loop
      total_mass_hot += BaryonField[DensNum][index];

    } // end temperature check
  } // end loop

  if( total_mass <= tiny_number){
    ENZO_FAIL("Grid_UpdateAverageAbundances: No mass in average abundances");
  } else{

    for(int im = 0; im <StellarYieldsNumberOfSpecies + 1; im++){
      if(mass_counter[im] <= 0.0  && mass_counter_hot[im] > 0.0){
        mass_counter[im] = mass_counter_hot[im];
      } else if (mass_counter[im] <= 0.0) {
        printf("Grid_UpdateAverageAbundances: Failing for species number %"ISYM" %"ESYM" %"ESYM" in list (where 0 = total metals)\n", im, mass_counter[im], mass_counter_hot[im]);
        ENZO_FAIL("Grid_UpdateAverageAbundances: No tracer species mass found on grid");
      }
    }
  }


  for(int im = 0; im < StellarYieldsNumberOfSpecies + 1; im ++){
    this->AveragedAbundances[im] = mass_counter[im] / total_mass;
    if(this->AveragedAbundances[im] <= 0.0){
        ENZO_FAIL("Averaged abundance ratio negative");
    }
  }

  delete [] temperature;
  return SUCCESS;
}

