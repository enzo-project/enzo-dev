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

int InitializeIMF(float *& data, const float & lower_mass, const float & upper_mass,
                  const int & IMFtype);
int InitializeDTD(float *& data);

float Ruiter_SNIa_DTD(float time, const int model);


unsigned_long_int mt_random(void);

int StarParticleIndividual_IMFInitialize(void){

  if (IMFData == NULL){
    InitializeIMF( IMFData, IndividualStarIMFLowerMassCutoff,
                            IndividualStarIMFUpperMassCutoff,
                            IndividualStarIMF);
  }

  if (IndividualStarPopIIIFormation && SecondaryIMFData == NULL){
    InitializeIMF( SecondaryIMFData, PopIIILowerMassCutoff, PopIIIUpperMassCutoff,
                                     3); // 3 == PopIII IMF
  }

  if (IndividualStarSNIaModel == 2){
    InitializeDTD( EventDTD );
  }

  return SUCCESS;
}

int InitializeDTD(float *& data){
 /* Initialize cumulative probability distribution from a
    delay type distribution */

  data = new float[IMF_TABLE_ENTRIES];

  const double min_time = log10(1.0E4); // yr
  const double max_time = log10(14.0E9); // yr
  const double dt = (max_time-min_time)/(double(IMF_TABLE_ENTRIES)-1);

  double t=0.0 ,t_prev = 0.0;
  const double inv_six = 1.0/6.0;
  double f_a = Ruiter_SNIa_DTD(POW(10.0,(min_time)), 0); // 0 = total
  data[0] = f_a;
  for (int i = 1; i < IMF_TABLE_ENTRIES; i++){
    t = POW(10.0, min_time + dt*i); // time in yr

    double f_b  = Ruiter_SNIa_DTD(t, 0);
    double f_ab = Ruiter_SNIa_DTD( 0.5*(t + t_prev), 0);
    data[i] = inv_six*(t - t_prev)*(f_a + 4.0*f_ab + f_b);

    data[i] *= IndividualStarSNIaFraction; // normalize
    f_a = f_b;
    t_prev = t;
  }

  // get cumulative dist
  for (int i = 1; i < IMF_TABLE_ENTRIES;i++){
    data[i] += data[i-1];
  }



  return SUCCESS;
}


int InitializeIMF(float *& data, const float & lower_mass, const float & upper_mass,
                                const int & IMFtype)
{

  data = new float[IMF_TABLE_ENTRIES];

  int i;
  float m, m0, dm, total_fn;

  dm = log10(upper_mass / lower_mass)/
       ((float) (IMF_TABLE_ENTRIES-1));
  m0 = log10(lower_mass);

  total_fn = 0; // will hold cumulative probability density function

  // use global parameters from global_data.h
  if (IMFtype == 0){ // Salpeter
    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);
      total_fn += POW(m, IndividualStarSalpeterSlope);
      data[i] = total_fn;
    } // end tabulate

  } else if (IMFtype == 1){ // Chabrier 2003

    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);
      total_fn += 0.158 * exp( - 0.5 * POW(log10(m)-log10(0.08),2.0) /
                                       POW(0.69,2) );
      data[i] = total_fn;

    } // end tabulate

  } else if (IMFtype == 2){ // Kroupa 2001
    for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
      m = POW(10.0, m0 + i*dm);

      if (m < 0.08){ // NEED TO CHECK UNITS !!!!
        total_fn += POW(m, IndividualStarKroupaAlpha1);
      } else if (m < 0.5){
        total_fn += POW(m, IndividualStarKroupaAlpha2);
      } else{
        total_fn += POW(m, IndividualStarKroupaAlpha3);
      }

      data[i] = total_fn;
    } // end tabulate

  } else if (IMFtype == 3){ // PopIII IMF
    const float CutoffExponent = 1.6;

    for (i = 0; i < IMF_TABLE_ENTRIES; i++){
      m = POW(10.0, m0 + i*dm);
      total_fn += POW(m, PopIIIInitialMassFunctionSlope) *
         exp(-POW((PopIIIStarMass/m), CutoffExponent));
      data[i] = total_fn;
    }
  } else{
    fprintf(stdout, "Error in StarParticleIndividual_IMFInitialize. Invalid IMF choice.\n");
    return FAIL;
  } // end IMF type check


  // Normalize cpdf to 1
  for (i = 0; i < IMF_TABLE_ENTRIES; i ++){
    data[i] /= data[IMF_TABLE_ENTRIES-1];
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
