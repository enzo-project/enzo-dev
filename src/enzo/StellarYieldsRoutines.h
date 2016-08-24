#ifndef __STELLAR_YIELD_ROUTINES_H
#define __STELLAR_YIELD_ROUTINES_H

#include "typedefs.h"

#ifdef DEFINE_STORAGE
# define ISEXTERN
#else
# define ISEXTERN extern
#endif


float StellarYields_SNIaYieldsByNumber(const int &atomic_number);

float StellarYieldsInterpolateYield(int yield_type, const float &M,
                                    const float &metallicity, int atomic_number);

int GetYieldIndex(const int &number_of_yields, const int &Z);

#endif