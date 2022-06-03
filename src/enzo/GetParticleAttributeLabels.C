/***********************************************************************
 *
 *  file    : GenerateParticleAttributeLabels.C
 *
 *  Author  : Andrew Emerick
 *  Date    : April 16, 2020
 *
 *  PURPOSE : Generate labels for particle attributes in a single,
 *            function since this gets used at multiple points in the code
 *
 *            An alternative to this is to construct the particle attribute
 *            labels at problem initialization the same way the baryon
 *            fields are built.
 *
 *            It was only after my 10th time adding something to this code
 *            that I got fed up doing it in 20 different places....
 *
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"

char* ChemicalSpeciesParticleLabel(const int &atomic_number);
char* IndividualStarTableIDLabel(const int &num);


void GetParticleAttributeLabels(std::vector<std::string> & ParticleAttributeLabel){

  if (NumberOfParticleAttributes == 0) return;

#ifdef WINDS
  const char *temp_labels[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x",
      "particle_jet_y", "particle_jet_z", "typeia_fraction"};
  for(int i = 0; i < 7; i++){
    ParticleAttributeLabel[i], temp_labels[i]);
  }
#else


  /* Assign labels conditionally based on runtime parameters */

  ParticleAttributeLabel[0] = "creation_time" ;

  ParticleAttributeLabel[1] = "dynamical_time";
  ParticleAttributeLabel[2] = "metallicity_fraction";


  if(STARMAKE_METHOD(INDIVIDUAL_STAR)){
    ParticleAttributeLabel[3] = "birth_mass";

    if(MultiMetals == 2 && !IndividualStarOutputChemicalTags){
      int ii = 0;
      for(int j = 0; j < StellarYieldsNumberOfSpecies; j++){
        ParticleAttributeLabel[4 + ii++] = ChemicalSpeciesParticleLabel(StellarYieldsAtomicNumbers[j]);
      }
      if (IndividualStarTrackAGBMetalDensity) ParticleAttributeLabel[4 + ii++] = "agb_metal_fraction";
      if (IndividualStarPopIIIFormation){
        ParticleAttributeLabel[4 + ii++] = "popIII_metal_fraction";
        ParticleAttributeLabel[4 + ii++] = "popIII_pisne_metal_fraction";

        if (IndividualStarPopIIISeparateYields){
          for(int j = 2; j < StellarYieldsNumberOfSpecies; j++){ // skip H and He
            ParticleAttributeLabel[4 + ii] = "popIII_";
            ParticleAttributeLabel[4 + ii++] += ChemicalSpeciesParticleLabel(StellarYieldsAtomicNumbers[j]);
          }
        }
      }

      if (IndividualStarTrackWindDensity){
        ParticleAttributeLabel[4 + ii++] = "intermediate_wind_metal_fraction";
        ParticleAttributeLabel[4 + ii++] = "massive_wind_metal_fraction";
      }

      if (IndividualStarTrackSNMetalDensity){
        ParticleAttributeLabel[4 + ii++] = "snia_metal_fraction";

        if (IndividualStarSNIaModel == 2){
          ParticleAttributeLabel[4 + ii++] = "snia_sch_metal_fraction";
          ParticleAttributeLabel[4 + ii++] = "snia_sds_metal_fraction";
          ParticleAttributeLabel[4 + ii++] = "snia_hers_metal_fraction";
        }

        ParticleAttributeLabel[4 + ii++] = "snii_metal_fraction";
      }
      if (IndividualStarRProcessModel) ParticleAttributeLabel[4 + ii++] = "rprocess_metal_fraction";


    } // endif multimetals

    if (IndividualStarSaveTablePositions){
      for(int ii = ParticleAttributeTableStartIndex; ii < NumberOfParticleAttributes; ii++){
        ParticleAttributeLabel[ii] = IndividualStarTableIDLabel(ii - ParticleAttributeTableStartIndex);
      }
    }
    ParticleAttributeLabel[NumberOfParticleAttributes-2] = "wind_mass_ejected";
    ParticleAttributeLabel[NumberOfParticleAttributes-1] = "sn_mass_ejected";

  } else { // not using individual star model

    if (StarMakerTypeIaSNe){
      ParticleAttributeLabel[3] = "typeia_fraction";
    }

  }

#endif

  return ;
}
