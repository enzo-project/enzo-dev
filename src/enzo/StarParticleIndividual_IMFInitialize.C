/**********************************************************
/
/ IMF LOOKUP TABLE INITIALIZATION FOR INDIVIDUAL STAR FORMATION
/
/
/ copied from: John Wise - StarParticlePopIII_IMFInitialize.C
/              in state as of Jan 2016
/              Copied to new function so does not break PopIII code
/ written by: Andrew Emerick 
/ date      : January, 2016
/ modified1 :
**********************************************************/

#include <stdlib.h>
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

void mt_init(unsigned_int seed);
//void random_init(unsigned_int seed); need to find this function (2/8/16 AJE)

unsigned_long_int mt_random(void);

int StarParticleIndividual_IMFInitialize(void)
{

  if (IMFData != NULL){
    return SUCCESS;
  }

  IMFData = new float[IMF_TABLE_ENTRIES];

  int i;
  float m, m0, dm, total_fn;

  dm = log10(IndividualStarIMFUpperMassCutoff / IndividualStarIMFMassFloor)/
       ((float) (IMF_TABLE_ENTRIES-1));
  m0 = log10(IndividualStarIMFMassFloor);

  total_fn = 0; // will hold cumulative probability density function

  // use global parameters from global_data.h
  if (IndividualStarIMF == 0){ // Salpeter
    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);
      total_fn += POW(m, IndividualStarSalpeterSlope);
      IMFData[i] = total_fn;
    } // end tabulate

  } else if (IndividualStarIMF == 1){ // Chabrier 2003

    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);
      total_fn += 0.158 * exp( - 0.5 * POW(log10(m)-log10(0.08),2.0) /
                                       POW(0.69,2) );
      IMFData[i] = total_fn;

    } // end tabulate
 
  } else if (IndividualStarIMF == 2){ // Kroupa 2001
    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);

      if (m < 0.08){ // NEED TO CHECK UNITS !!!!
        total_fn += POW(m, IndividualStarKroupaAlpha1);
      } else if (m < 0.5){
        total_fn += POW(m, IndividualStarKroupaAlpha2);
      } else{
        total_fn += POW(m, IndividualStarKroupaAlpha3);
      }

      IMFData[i] = total_fn;
    } // end tabulate

  } else{
    fprintf(stdout, "Error in StarParticleIndividual_IMFInitialize. Invalid IMF choice.\n");
    return FAIL;
  } // end IMF type check


  // Normalize cpdf to 1
  for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
    IMFData[i] /= IMFData[IMF_TABLE_ENTRIES-1];
  }

  // Initialize random number generator. If restart, call it the
  // number of times + 1 to return to state before restart.
  // As of June 2016, we are breaking repeatability to have each
  // processor with their own random number generator (otherwise
  // SF regions will be identical across the simulation)

  if (IndividualStarIMFSeed == INT_UNDEFINED){
    mt_init( time(NULL) + MyProcessorNumber );
  }
  else{
    mt_init(IndividualStarIMFSeed + MyProcessorNumber);
  }

  unsigned_long_int trash;
  for (i = 0; i < 1 + IndividualStarIMFCalls; i++){
    trash = mt_random();
  }

  return SUCCESS;
}


