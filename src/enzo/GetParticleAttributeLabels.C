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

void GetParticleAttributeLabels(char * ParticleAttributeLabel){


#ifdef WINDS
  const char *temp_labels[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x",
      "particle_jet_y", "particle_jet_z", "typeia_fraction"};
  for(int i = 0; i < 7; i++){
    strcpy(&ParticleAttributeLabel[i], temp_labels[i]);
  }
#else


  /* Assign labels conditionally based on runtime parameters */

  strcpy(&ParticleAttributeLabel[0],"creation_time");
  strcpy(&ParticleAttributeLabel[1],"dynamical_time");
  strcpy(&ParticleAttributeLabel[2],"metallicity_fraction");

  if(STARMAKE_METHOD(INDIVIDUAL_STAR)){
    strcpy(&ParticleAttributeLabel[3],"birth_mass");

    if(MultiMetals == 2 && !IndividualStarOutputChemicalTags){
      int ii = 0;
      for(ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
        strcpy(&ParticleAttributeLabel[4 + ii],ChemicalSpeciesParticleLabel(StellarYieldsAtomicNumbers[ii]));
      }
      if (IndividualStarTrackAGBMetalDensity) strcpy(&ParticleAttributeLabel[4 + ii++],"agb_metal_fraction");
      if (IndividualStarPopIIIFormation){
        strcpy(&ParticleAttributeLabel[4 + ii++],"popIII_metal_fraction");
        strcpy(&ParticleAttributeLabel[4 + ii++],"popIII_pisne_metal_fraction");
      }

      if (IndividualStarTrackSNMetalDensity){
        strcpy(&ParticleAttributeLabel[4 + ii++],"snia_metal_fraction");

        if (IndividualStarSNIaModel == 2){
          strcpy(&ParticleAttributeLabel[4 + ii++], "snia_sch_metal_fraction");
          strcpy(&ParticleAttributeLabel[4 + ii++], "snia_sds_metal_fraction");
          strcpy(&ParticleAttributeLabel[4 + ii++], "snia_hers_metal_fraction");
        }

        strcpy(&ParticleAttributeLabel[4 + ii++],"snii_metal_fraction");
      }
      if (IndividualStarRProcessModel) strcpy(&ParticleAttributeLabel[4 + ii++],"rprocess_metal_fraction");


    } // endif multimetals

    if (IndividualStarSaveTablePositions){
      for(int ii = ParticleAttributeTableStartIndex; ii < NumberOfParticleAttributes; ii++){
        strcpy(&ParticleAttributeLabel[ii],IndividualStarTableIDLabel(ii - ParticleAttributeTableStartIndex));
      }
    }
    strcpy(&ParticleAttributeLabel[NumberOfParticleAttributes-2],"wind_mass_ejected");
    strcpy(&ParticleAttributeLabel[NumberOfParticleAttributes-1],"sn_mass_ejected");

  } else { // not using individual star model

    strcpy(&ParticleAttributeLabel[3],"typeia_fraction");
  }

#endif

  return ;
}
