/***********************************************************************
/
/  COMPUTE PHOTON PRODUCTION RATES FOR (NON)ROTATING POPIII STARS
/
/  written by: Danielle Skinner
/  date:       September, 2021
/  modified1:
/
/  ---------- SPECIES --------
/  0 : HI
/  1 : HeI
/  2 : HeII
/  
/ Model comes from Murphy et al. 2021. Only valid for stars with masses 
/ between 9 - 120 Msun. H2 not included. The fits here are up until the 
/ star is at 80% of its lifetime, after which the photon production rate
/ is assumed to be constant.
/ 
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "phys_constants.h"
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

int j, bin, nvalue;
float frac, a, b, c, Mass, age;
double Q, Q_lower, Q_upper;

const int nentries = 9;
const float mass_list[] = {9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 60.0, 85.0, 120.0}; 

int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);

int find_element(float element, const float array[], int n){
    nvalue = -1;
  for (j = 0; j < n; j++){
      if(array[j] == element){
          nvalue = j;
          break;
      }
  }
  return nvalue;
}

double InterpolateQ(float Mass, float age, float *parameters_a, float *parameters_b, float *parameters_c){
    float age2 = age*age;

    if (find_element(Mass, (float*)mass_list, nentries) != -1) {
        Q = (parameters_a[nvalue] * age2 + parameters_b[nvalue] * age + parameters_c[nvalue]);
        return Q;
      }
    else {
        bin = search_lower_bound((float*)mass_list, Mass, 0, nentries, nentries);
        frac = (Mass - mass_list[bin]) / (mass_list[bin+1] - mass_list[bin]);
        Q_lower = ( (parameters_a[bin] * age2) + (parameters_b[bin] * age) + parameters_c[bin] );
        Q_upper = ( (parameters_a[bin+1] * age2) + (parameters_b[bin+1] * age) + parameters_c[bin+1] );

        Q = Q_lower + frac * (Q_upper - Q_lower);
        return Q;
      }
}

double CalculateRotationalPhotonRates(float Mass, float age, int species){

    if (species == 0) {
        static float parameters_a[] = {-0.3014, 0.03601, 0.1357, 0.2404, 0.2116, 0.1807, 0.1118, 0.0901, 0.0683};
        static float parameters_b[] = {0.66228, 0.40061, 0.26264, 0.15633, 0.17335, 0.177001, 0.21019, 0.19481, 0.19028};
        static float parameters_c[] = {47.42321, 47.82289, 48.15684, 48.52598, 48.96249, 49.23924, 49.58557, 49.85428, 50.09629};

        Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);
      
    }

    else if (species == 1) {
        static float parameters_a[] = {-0.9481, -0.6126, -0.3906, -0.1218, -0.1451, -0.1122, -0.4014, -0.4259, -0.4637};
        static float parameters_b[] = {1.19688, 0.99590, 0.65148, 0.32491, 0.28988, 0.23664, 0.37460, 0.33801, 0.31139};
        static float parameters_c[] = {46.48224, 46.98303, 47.47341,  47.98433, 48.52012, 48.84545, 49.22893, 49.52746, 49.79099};
    
        Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);

    }

    else if (species == 2) {
        static float parameters_a[] = {-2.8125, -2.5140, -1.9560, -1.2121, -1.2286, -1.0128, -1.9777, -2.0174, -2.0865};
        static float parameters_b[] = {2.73813, 2.74062, 1.80823, 0.83063, 0.63568, 0.41093, 0.86035, 0.75522, 0.63998};
        static float parameters_c[] = {43.43041, 44.22116, 45.17038, 46.10727, 46.95021, 47.42831, 47.93338, 48.32898, 48.66499};
    
        Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);

    }
    return Q;
}

double CalculateNonRotationalPhotonRates(float Mass, float age, int species){

    if (species == 0) {
        static float parameters_a[] = {-0.35938, -0.03421, 0.05179, 0.17966, 0.18772, 0.17035, 0.15019, 0.12873, 0.12268};
        static float parameters_b[] = {0.64365, 0.38223, 0.26300, 0.13555, 0.13217, 0.14523, 0.14671, 0.14156, 0.12346};
        static float parameters_c[] = {47.444944, 47.830253, 48.174224, 48.544643, 48.979915, 49.254436, 49.601645, 49.868288, 50.109797};
    
      Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);

    }

    else if (species == 1) {
        static float parameters_a[] = {-1.05193, -0.53651, -0.46799, -0.15949, -0.08082, -0.11586, -0.17175, -0.22875, -0.27835};
        static float parameters_b[] = {1.19633, 0.95788, 0.65528, 0.30155, 0.21199, 0.22530, 0.23914, 0.24920, 0.25542};
        static float parameters_c[] = {46.522791, 46.981776, 47.505882, 48.018253, 48.553852, 48.872077, 49.263237, 49.555736, 49.814517};
    
      Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);

    }

    else if (species == 2) {
        static float parameters_a[] = {-3.05209, -2.01521, -2.01975, -1.18593, -0.90444, -1.00130, -1.17825, -1.35364, -1.54739};
        static float parameters_b[] = {2.79275, 2.65228, 1.82661, 0.80239, 0.45104, 0.46494, 0.51727, 0.57479, 0.65844};
        static float parameters_c[] = {43.524763, 44.193070, 45.247383, 46.187693, 47.034429, 47.491437, 48.024623, 48.402114, 48.718482};
    
      Q = InterpolateQ(Mass, age, parameters_a, parameters_b, parameters_c);

    }

  return Q;

}