#ifndef __STELLAR_YIELD_ROUTINES_H
#define __STELLAR_YIELD_ROUTINES_H

// #include <hdf5.h>


//#include "typedefs.h"

#ifdef DEFINE_STORAGE
# define ISEXTERN
#else
# define ISEXTERN extern
#endif


float StellarYields_ScaledSolarMassFractionByNumber(const float &metallicity,
                                                    const int   &atomic_number,
                                                    const int table = 0);


float StellarYields_SolarAbundancesByNumber(const int &atomic_number);
float StellarYields_SolarAbundancesByNumber(const int &atomic_number, const int table);
float StellarYields_SolarAbundancesByNumber_Asplund(const int &atomic_number);
float StellarYields_SolarAbundancesByNumber_Lodders(const int &atomic_number);


float StellarYields_MMW(const int &atomic_number);

float StellarYields_SNIaYieldsByNumber(const int &atomic_number);

// float StellarYields_PopIIIYieldsByNumber(const int &atomic_number);

float StellarYieldsInterpolatePopIIIYield(const int &i, const float &M, int atomic_number);

float StellarYieldsInterpolateYield(int yield_type, const float &M,
                                    const float &metallicity, int atomic_number);

int StellarYieldsGetYieldTablePosition(int &i, int &j,
                                        const float &M, const float &metallciity);

int StellarYieldsGetYieldTablePosition(const StellarYieldsDataType & table,
                                       int &i, int &j, const float &M, const float &metallicity,
                                       int extrapolate_yields_switch = 0);

float StellarYieldsInterpolateYield(int yield_type,
                                    const int &i, const int &j,
                                    const float &M, const float &metallicity, int atomic_number);

int StellarYieldsGetPopIIIYieldTablePosition(int &i, const float &M);

int GetYieldIndex(const int &number_of_yields, const int &Z);


#ifdef NEWYIELDTABLES
int read_dataset(hid_t file_id, const char *dset_name, double *buffer);
void initialize_table_1D(StellarYieldsDataType table);
void fill_table(StellarYieldsDataType *table,
                   std::string filename,
                   std::string dname);
#endif


#endif
