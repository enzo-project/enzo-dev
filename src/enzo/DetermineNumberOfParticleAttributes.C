/***********************************************************************
 *
 *  file    : DetermineNumberOfParticleAttributes.C
 *
 *  Author  : Andrew Emerick
 *  Date    : May 4, 2020
 *
 *  PURPOSE : Go through all of the parameter settings to determine the number
 *            of particle attributes. Just returns this number rather than
 *            setting the globa value itself.
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


int DetermineNumberOfAbundanceAttributes(void){
  /// an abridged version of DetermineNumberOfParticleAttributes
  // to just get the number of abundances

  int n = StellarYieldsNumberOfSpecies;

  if (IndividualStarTrackAGBMetalDensity) n++;
  if (IndividualStarPopIIIFormation) n += 2;
  if (IndividualStarPopIIISeparateYields) n += (StellarYieldsNumberOfSpecies -2);
  if (IndividualStarTrackWindDensity) n+= 1;
  if (IndividualStarTrackSNMetalDensity){
    n += 2;
    if (IndividualStarSNIaModel == 2) n += 3;
  }
  if (IndividualStarRProcessModel) n += 1;

  return n;
}


int DetermineNumberOfParticleAttributes(void){

  int n=0;

  if (StarParticleCreation || StarParticleFeedback) {


    // first three attributes are:
    //   1) Formation time
    //   2) Dynamical time / Life time
    //   3) Metallicity

    n = 3;

    if (StarMakerTypeIaSNe) n++;
    if (StarMakerTypeIISNeMetalField) n++;

    if (STARMAKE_METHOD(INDIVIDUAL_STAR)){
      n++; // birth mass

      if (MultiMetals){

        if(MultiMetals == 2 && !IndividualStarOutputChemicalTags){
          n += DetermineNumberOfAbundanceAttributes();
        }


      } // end multi metals
      if (IndividualStarSaveTablePositions){
        ParticleAttributeTableStartIndex = n;
        n += NumberOfParticleTableIDs;
      }
      n += 2; // counters for mass loss
    }

  } else {
    n = 0;
  }

  return n;
}
