#ifndef __STELLAR_YIELD_ROUTINES_H
#define __STELLAR_YIELD_ROUTINES_H

#include "typedefs.h"

#ifdef DEFINE_STORAGE
# define ISEXTERN
#else
# define ISEXTERN extern
#endif


float StellarYields_SNIaYieldsByNumber(const int &atomic_number);

// float StellarYields_PopIIIYieldsByNumber(const int &atomic_number);

float StellarYieldsInterpolatePopIIIYield(const int &i, const float &M, int atomic_number);

float StellarYieldsInterpolateYield(int yield_type, const float &M,
                                    const float &metallicity, int atomic_number);

int StellarYieldsGetYieldTablePosition(int &i, int &j,
                                        const float &M, const float &metallciity);

float StellarYieldsInterpolateYield(int yield_type,
                                    const int &i, const int &j,
                                    const float &M, const float &metallicity, int atomic_number);

int GetYieldIndex(const int &number_of_yields, const int &Z);

#endif
