/*************************************************************************
/
/ STELLAR YIELDS DATA: CHEMICAL YIELDS, MECHANICAL LUMINOSITY, MASS EJECTA
/
/ Written by: A. Emerick
/ Date      : 3/6/16
/
/   WARNING: In every case where interpolation is over metallicity, we institute
/            a 'metallicity floor' whereby the star's metallicity is set to Z_min
/            during interpolation for that given table, whenever Z_star < Z_min
/            for that given table. This *may* lead to potential inconsistencies if
/            Z_min is different in each table, but really that means certain
/            properties will be *better* than others in this regime. This is also
/            applied to the stellar properties tables as well.
/
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "StellarYieldsRoutines.h"
#include "StarParticleData.h"

int search_lower_bound(float *arr, float value, int low, int high,
                       int total);


float StellarYields_SolarAbundancesByNumber(const int &atomic_number){
  /* For backwards compatability */
  return StellarYields_SolarAbundancesByNumber(atomic_number, 0);
}

float StellarYields_SolarAbundancesByNumber(const int &atomic_number,
                                            const int table){

    float value = 0.0;

    if (table == 0){ // use Asplund 2009
        value = StellarYields_SolarAbundancesByNumber_Asplund(atomic_number);
    } else if (table == 1){
        value = StellarYields_SolarAbundancesByNumber_Lodders(atomic_number);
    } else {
      ENZO_VFAIL("Error loading solar abundance table - incorrect table number\n.");
    }

    return value;
}

float StellarYields_SolarAbundancesByNumber_Lodders(const int &atomic_number){
   /* From Lodders+2003 */

  float abund = 0.0;

  switch(atomic_number){
    case -2: abund = 0.2741; break; // solar He mass fraction
    case -1: abund = 0.7110; break; // solar H mass fraction
    case  0: abund = 0.0149; break; // Solar Mass fraction of metals - NOT ABUNDANCE
    case  1: abund = 12.00; break;    case  2: abund = 10.984; break;
    case  3: abund =  3.35; break;    case  4: abund =  1.48; break;
    case  5: abund =  2.85; break;    case  6: abund =  8.46; break;
    case  7: abund =  7.90; break;    case  8: abund =  8.76; break;
    case  9: abund =  4.53; break;    case 10: abund =  7.95; break;
    case 11: abund =  6.37; break;    case 12: abund =  7.62; break;
    case 13: abund =  6.54; break;    case 14: abund =  7.61; break;
    case 15: abund =  5.54; break;    case 16: abund =  7.26; break;
    case 17: abund =  5.33; break;    case 18: abund =  6.62; break;
    case 19: abund =  5.18; break;    case 20: abund =  6.41; break;
    case 21: abund =  3.15; break;    case 22: abund =  5.00; break;
    case 23: abund =  4.07; break;    case 24: abund =  5.72; break;
    case 25: abund =  5.58; break;    case 26: abund =  7.54; break;
    case 27: abund =  4.98; break;    case 28: abund =  6.29; break;
    case 29: abund =  4.34; break;    case 30: abund =  4.70; break;
    case 31: abund =  3.17; break;    case 32: abund =  3.70; break;
    case 33: abund =  2.40; break;    case 34: abund =  3.43; break;
    case 35: abund =  2.67; break;    case 36: abund =  3.36; break;
    case 37: abund =  2.43; break;    case 38: abund =  2.99; break;
    case 39: abund =  2.28; break;    case 40: abund =  2.67; break;
    case 41: abund =  1.49; break;    case 42: abund =  2.03; break;
    // case 43: abund = -?; break; Tc - go to default
    case 44: abund =  1.89; break;    case 45: abund =  1.18; break;
    case 46: abund =  1.77; break;    case 47: abund =  1.30; break;
    case 48: abund =  1.81; break;    case 49: abund =  0.87; break;
    case 50: abund =  2.19; break;    case 51: abund =  1.14; break;
    case 52: abund =  2.30; break;    case 53: abund =  1.61; break;
    case 54: abund =  2.35; break;    case 55: abund =  1.18; break;
    case 56: abund =  2.25; break;    case 57: abund =  1.25; break;
    case 58: abund =  1.68; break;    case 59: abund =  0.85; break;
    case 60: abund =  1.54; break;    case 61: abund =  0.00; break;
    case 62: abund =  1.02; break;    case 63: abund =  0.60; break;
    case 64: abund =  1.13; break;    case 65: abund =  0.38; break;
    case 66: abund =  1.21; break;    case 67: abund =  0.56; break;
    case 68: abund =  1.02; break;    case 69: abund =  0.18; break;
    case 70: abund =  1.01; break;    case 71: abund =  0.16; break;
    case 72: abund =  0.84; break;    case 73: abund = -0.06; break;
    case 74: abund =  0.72; break;    case 75: abund =  0.33; break;
    case 76: abund =  1.44; break;    case 77: abund =  1.42; break;
    case 78: abund =  1.75; break;    case 79: abund =  0.91; break;
    case 80: abund =  1.23; break;    case 81: abund =  0.88; break;
    case 82: abund =  2.13; break;    case 83: abund =  0.76; break;

    default:
      ENZO_FAIL("Failure in StellarYields_SolarAbundancesByNumber_Lodders: Wrong atomic number");
  }

  return abund;
}



float StellarYields_SolarAbundancesByNumber_Asplund(const int &atomic_number){
 /* ----------------------------------------------------------------
  * StellarYields_SolarAbundancesByNumber
  * ----------------------------------------------------------------
  * Given atomic number of an element, returns the solar abundances
  * of that element as reported in Table 1 of
  * Asplund et. al. 2009 (ARAA 47:481-522). All abundances are
  * the recorded photosphere abundances when available, but meteoric
  * abundances are used instead when these are unavailable (for
  * As, Se, Br, Cd, Sb, Te, I, Cs, Ta, Re, Pt, Hg, Bi). Abundances
  * here are reported as-is in table, of the form:
  *
  * log(e_x) = log(N_x / N_H) + 12.0
  *
  * Where the LHS is the value reported in the table for element X,
  * N_x is the abundance of element X, and N_H is abundance of
  * Hydrogen.
  *
  * Using 0 returns the metal mass fraction
  * Using -1 returns the H mass fraction
  * Using -2 returns the He mass fraction
  * ---------------------------------------------------------------*/

  float abund = 0.0;

  switch(atomic_number){
    case -2: abund = 1.0 - 0.7381 - 0.0134; break; // He fraction by mass
    case -1: abund = 0.7381; break; // H fraction by mass
    case  0: abund = 0.0134; break; // Solar Mass fraction of metals - NOT ABUNDANCE
    case  1: abund = 12.00; break;    case  2: abund = 10.93; break;
    case  3: abund =  1.05; break;    case  4: abund =  1.38; break;
    case  5: abund =  2.70; break;    case  6: abund =  8.43; break;
    case  7: abund =  7.83; break;    case  8: abund =  8.69; break;
    case  9: abund =  4.56; break;    case 10: abund =  7.93; break;
    case 11: abund =  6.24; break;    case 12: abund =  7.60; break;
    case 13: abund =  6.45; break;    case 14: abund =  7.51; break;
    case 15: abund =  5.41; break;    case 16: abund =  7.12; break;
    case 17: abund =  5.50; break;    case 18: abund =  6.40; break;
    case 19: abund =  5.03; break;    case 20: abund =  6.34; break;
    case 21: abund =  3.15; break;    case 22: abund =  4.95; break;
    case 23: abund =  3.93; break;    case 24: abund =  5.64; break;
    case 25: abund =  5.43; break;    case 26: abund =  7.50; break;
    case 27: abund =  4.99; break;    case 28: abund =  6.22; break;
    case 29: abund =  4.19; break;    case 30: abund =  4.56; break;
    case 31: abund =  3.04; break;    case 32: abund =  3.65; break;
    case 33: abund =  2.30; break;    case 34: abund =  3.34; break;
    case 35: abund =  2.54; break;    case 36: abund =  3.25; break;
    case 37: abund =  2.52; break;    case 38: abund =  2.87; break;
    case 39: abund =  2.21; break;    case 40: abund =  2.58; break;
    case 41: abund =  1.46; break;    case 42: abund =  1.88; break;
    // case 43: abund =  0.00; break; // Tc is not natural
    case 44: abund =  1.75; break;    case 45: abund =  0.91; break;
    case 46: abund =  1.57; break;    case 47: abund =  0.94; break;
    case 48: abund =  1.71; break;    case 49: abund =  0.80; break;
    case 50: abund =  2.04; break;    case 51: abund =  1.01; break;
    case 52: abund =  2.18; break;    case 53: abund =  1.55; break;
    case 54: abund =  2.24; break;    case 55: abund =  1.08; break;
    case 56: abund =  2.18; break;    case 57: abund =  1.10; break;
    case 58: abund =  1.58; break;    case 59: abund =  0.72; break;
    case 60: abund =  1.42; break;    case 61: abund =  0.00; break;
    case 62: abund =  0.96; break;    case 63: abund =  0.52; break;
    case 64: abund =  1.07; break;    case 65: abund =  0.30; break;
    case 66: abund =  1.10; break;    case 67: abund =  0.48; break;
    case 68: abund =  0.92; break;    case 69: abund =  0.10; break;
    case 70: abund =  0.84; break;    case 71: abund =  0.10; break;
    case 72: abund =  0.85; break;    case 73: abund = -0.12; break;
    case 74: abund =  0.85; break;    case 75: abund =  0.26; break;
    case 76: abund =  1.40; break;    case 77: abund =  1.38; break;
    case 78: abund =  1.62; break;    case 79: abund =  0.92; break;
    case 80: abund =  1.17; break;    case 81: abund =  0.90; break;
    case 82: abund =  1.75; break;    case 83: abund =  0.65; break;

    default:
      ENZO_FAIL("Failure in StellarYields_SolarAbundancesByNumber_Asplund: Wrong atomic number");
  }

  return abund;
}

float StellarYields_MMW(const int &atomic_number){

  float A = 0.0;

  switch(atomic_number){
    case 1: A = 1.007900; break;     case 2: A = 4.002600; break;
    case 3: A = 6.941000; break;     case 4: A = 9.012200; break;
    case 5: A = 10.811000; break;     case 6: A = 12.010700; break;
    case 7: A = 14.006700; break;     case 8: A = 15.999400; break;
    case 9: A = 18.998400; break;     case 10: A = 20.179700; break;
    case 11: A = 22.989700; break;     case 12: A = 24.305000; break;
    case 13: A = 26.981500; break;     case 14: A = 28.085500; break;
    case 15: A = 30.973800; break;     case 16: A = 32.065000; break;
    case 17: A = 35.453000; break;     case 18: A = 39.948000; break;
    case 19: A = 39.098000; break;     case 20: A = 40.078000; break;
    case 21: A = 44.955912; break;     case 22: A = 47.867000; break;
    case 23: A = 50.941500; break;     case 24: A = 51.996100; break;
    case 25: A = 54.938045; break;     case 26: A = 55.845000; break;
    case 27: A = 58.933195; break;     case 28: A = 58.693400; break;
    case 29: A = 63.546000; break;     case 30: A = 65.380000; break;
    case 31: A = 69.723000; break;     case 32: A = 72.640000; break;
    case 33: A = 74.921600; break;     case 34: A = 78.960000; break;
    case 35: A = 79.904000; break;     case 36: A = 83.798000; break;
    case 37: A = 85.467800; break;     case 38: A = 87.620000; break;
    case 39: A = 88.905850; break;     case 40: A = 91.224000; break;
    case 41: A = 92.906380; break;     case 42: A = 95.960000; break;
    case 43: A = 97.907200; break;     case 44: A = 101.070000; break;
    case 45: A = 102.905500; break;     case 46: A = 106.420000; break;
    case 47: A = 107.868200; break;     case 48: A = 112.411000; break;
    case 49: A = 114.818000; break;     case 50: A = 118.710000; break;
    case 51: A = 121.760000; break;     case 52: A = 127.600000; break;
    case 53: A = 126.904470; break;     case 54: A = 131.293000; break;
    case 55: A = 132.9054519; break;     case 56: A = 137.327000; break;
    case 57: A = 138.905470; break;     case 58: A = 140.116000; break;
    case 59: A = 140.907650; break;     case 60: A = 144.242000; break;
    case 61: A = 145.000000; break;     case 62: A = 150.360000; break;
    case 63: A = 151.964000; break;     case 64: A = 157.250000; break;
    case 65: A = 158.925350; break;     case 66: A = 162.500000; break;
    case 67: A = 164.930320; break;     case 68: A = 167.259000; break;
    case 69: A = 168.934210; break;     case 70: A = 173.054000; break;
    case 71: A = 174.966800; break;     case 72: A = 178.490000; break;
    case 73: A = 180.947880; break;     case 74: A = 183.840000; break;
    case 75: A = 186.207000; break;     case 76: A = 190.230000; break;
    case 77: A = 192.217000; break;     case 78: A = 195.084000; break;
    case 79: A = 196.966569; break;     case 80: A = 200.590000; break;
    case 81: A = 204.383300; break;     case 82: A = 207.200000; break;

    default:
      ENZO_FAIL("Failure in StellarYields_AtomicMassByNumber: Wrong atomic number");
  }

  return A;
}

float StellarYields_AtomicMassByNumber(const int &atomic_number){
  return StellarYields_MMW(atomic_number) * AMU_CGS;
}

float StellarYields_ScaledSolarMassFractionByNumber(const float &metallicity,
                                                    const int   &atomic_number,
                                                    const int table // default 0
                                                   ){

  // Constants below taken from values reported in
  // is the source for the solar abundances used to scale
  // [see section 3.12 of this work]

  const float solar_H_mass_fraction = StellarYields_SolarAbundancesByNumber(-1, table);
  const float solar_metallicity     = StellarYields_SolarAbundancesByNumber(0, table);

  float Z = metallicity / solar_metallicity;

  if (atomic_number <= 2){
    ENZO_FAIL("Failure in StellarYields_ScaledSolarMassFractionByNumber: cannot set metallicity scaled abundances for H or He");
  }

  float solar_abundance;
  solar_abundance = StellarYields_SolarAbundancesByNumber(atomic_number, table);

  // Asplund abundances are reported as log(e_x) = log(N_x/N_H) + 12.0
  // where LHS is the value in table, need to remove 12 scaling to get
  // actual abundance
  float e_x = POW(10.0, solar_abundance - StellarYields_SolarAbundancesByNumber(1, table));

  // solar mass fraction
  float f_x = e_x * solar_H_mass_fraction * (StellarYields_MMW(atomic_number) /
                                             StellarYields_MMW(1));

  // scaled solar mass fraction
  return f_x * Z;

//  return solar_H_mass_fraction * POW(10.0, solar_abundance - 12.0) * Z;
}

/*
float StellarYields_PopIIIYieldsByNumber(const int &atomic_number){
// Adopted as Z = 0 yields from Nomoto et. al. 2006


  float yield = -1.0;

  if (atomic_number < 0){
    return 1.0; // total mass as summed from table
  }

  if (atomic_number == 0){
    return 1.0; // total mass in metals
  }

  switch(atomic_number){
    case  1 : yield = ; break;
    case  2 : yield = ; break;
    default:
      yield = 1.0;
  }

  return yield;
}
*/

float StellarYields_SNIaYieldsByNumber(const int &atomic_number){
  /* -------------------------------------------------------------
   * StellarYields_SNIaYieldsByNumber
   * -------------------------------------------------------------
   * Given the atomic number of a desired element, returns the SNIa
   * yields for a given Type Ia supernova as given in the
   * Thieleman et. al. 1986 table (see Table 5)
   *
   * Table 5 in Thielemann et. al. 1986 has yields for multiple
   * isotopes for each species.
   * -------------------------------------------------------------*/

  float yield = -1.0;

  if (atomic_number <= 0) {
    return 1.2447714757 ; // total mass as summed from table (all mass is in metals)
  }


  if (atomic_number <= 5 || atomic_number >= 33 ){
    return 0.0; // nothing produced here (not included in table)
  }

  // yields are presented below as sum over all isotpes
  // in table in units of Msun
  switch(atomic_number){
    case  6 : yield = 5.0E-2 + 4.5E-13; break;
    case  7 : yield = 2.7E-9 + 4.4E-9 ; break;
    case  8 : yield = 1.3E-1 + 1.1E-10 + 1.7E-12; break;
    case  9 : yield = 2.5E-13; break;
    case 10 : yield = 1.8E-3 + 1.1E-8  + 2.5E-3; break;
    case 11 : yield = 1.8E-6; break;
    case 12 : yield = 1.6E-6 + 5.8E-6 + 4.0E-6; break;
    case 13 : yield = 4.4E-4; break;
    case 14 : yield = 1.5E-1 + 3.0E-4 + 3.4E-3; break;
    case 15 : yield = 1.4E-4; break;
    case 16 : yield = 8.2E-2 + 7.2E-4 + 1.5E-3 + 2.5E-8; break;
    case 17 : yield = 1.2E-4 + 2.8E-5; break;
    case 18 : yield = 1.7E-2 + 1.2E-3; break;
    case 19 : yield = 9.9E-5 + 6.6E-6; break;
    case 20 : yield = 1.5E-2 + 3.6E-5 + 4.2E-8 + 1.8E-5 + 1.3E-9 + 5.7E-12; break;
    case 21 : yield = 1.6E-7; break;
    case 22 : yield = 1.9E-5 + 3.1E-7 + 2.0E-4 + 9.3E-6 + 1.6E-6; break;
    case 23 : yield = 5.0E-9 + 2.8E-5; break;
    case 24 : yield = 2.3E-4 + 5.2E-3 + 6.6E-4 + 3.8E-5; break;
    case 25 : yield = 6.7E-3; break;
    case 26 : yield = 9.0E-2 + 6.3E-1 + 2.2E-2 + 2.5E-4; break;
    case 27 : yield = 7.3E-4; break;
    case 28 : yield = 1.3E-2 + 1.4E-2 + 2.4E-4 + 5.1E-3 + 2.6E-7; break;
    case 29 : yield = 2.0E-6 + 8.5E-6; break;
    case 30 : yield = 1.3E-5 + 1.9E-5 + 8.2E-8 + 3.5E-7 + 1.0E-9; break;
    case 31 : yield = 1.0E-7 + 6.1E-9; break;
    case 32 : yield = 8.4E-7 + 5.7E-8 + 8.1E-11 + 1.8E-8; break;

    default:
      ENZO_FAIL("Failure to find atomic number in SNIa yield computation");
  }

  return yield;
}


int StellarYieldsGetYieldTablePosition(const StellarYieldsDataType & table,
                                       int &i, int &j, const float &M, const float &metallicity,
                                       int extrapolate_yields_switch){

  int width, bin_number;

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table.M[0]) || (M > table.M[table.Nm - 1])){

    if( IndividualStarExtrapolateYields || extrapolate_yields_switch ){
      // use nearest mass bin
      const float fudge = 0.000001; // mass needs to be within bounds
      interp_M  = (M < table.M[0]) ? table.M[0]*(1.0+fudge): table.M[table.Nm - 1] * (1.0-fudge);

    } else{
      printf("StellarYieldsInterpolateYield: Mass out of bounds\n");
      printf("M = %"ESYM" for minimum M = %"ESYM" and maximum M = %"ESYM"\n", M, table.M[0], table.M[table.Nm-1]);
      return FAIL;
    }

  }

  if( (Z < table.Z[0]) || (Z > table.Z[table.Nz - 1])){

    if ( Z < table.Z[0] ){ // WARNING: see statement at top of file
      Z = table.Z[0];
    } else {

      printf("StellarYieldsInterpolateYield: Metallicity out of bounds\n");
      printf("Z = %"ESYM" for minimum Z = %"ESYM" and maximum Z = %"ESYM"\n", Z, table.Z[0], table.Z[table.Nz-1]);

      return FAIL;
    }
  }

  i = search_lower_bound(table.M, interp_M, 0, table.Nm, table.Nm);
  j = search_lower_bound(table.Z, Z, 0, table.Nz, table.Nz);

  if( (interp_M > table.M[i+1]) || interp_M < (table.M[i])){
    printf("interp_m = %"ESYM" i = %"ISYM" j = %"ISYM"\n", interp_M, i, j);
    ENZO_FAIL("FAIURE IN STELLAR YIELDS TABLE INTERPOLATION");
  }

  return SUCCESS;
}

int StellarYieldsGetYieldTablePosition(int &i, int &j,
                                        const float &M, const float &metallicity){
/* ------------------------------------------------------------------------
 * StellarYieldsGetYieldTablePosition
 * ------------------------------------------------------------------------
 * Interpolation function which finds the position in the yield table (i and j)
 *
 * This function tries to guess which yield table is relevant. For backwards
 * compatability, this has a few extra if statements
 * -------------------------------------------------------------------------*/
  /* interpolate table */
  StellarYieldsDataType * table;


#ifdef NEWYIELDTABLES
  /* In new yield tables */

  if (M <= IndividualStarAGBThreshold){

    if (StellarYieldsAGBData.M == NULL){
      //  This is done for backwards compatability
      //  In original tables, winds contained AGB stars
      table = &StellarYieldsWindData;

    } else {
      // otherwise we are doig the new style
      table = &StellarYieldsAGBData;
    }

  } else{

    if (StellarYieldsMassiveStarData.M != NULL){
      // this is the OLD style for backwards compatability
      if (M <= StellarYieldsSNData.M[StellarYieldsSNData.Nm-1]){
        table = &StellarYieldsSNData;
      } else if (( M >= StellarYieldsMassiveStarData.M[0]) &&
                 ( M <= StellarYieldsMassiveStarData.M[StellarYieldsMassiveStarData.Nm-1])){
         // use massive star data from PARSEC models (massive stars only! and wind only!)

         table = &StellarYieldsMassiveStarData;
       }

    } else {
      // in the new model, both the SN table and the winds tables have
      // the same mass bins. In this case, it doesn't matter which we pick
      table = &StellarYieldsWindData;
    }
  }

#else

  if (M <= StellarYieldsSNData.M[StellarYieldsSNData.Nm - 1]){
    // if mass is on NuGrid data set (M <= 25) use this table
    table = &StellarYieldsSNData;
  } else if (( M >= StellarYieldsMassiveStarData.M[0]) &&
            ( M <= StellarYieldsMassiveStarData.M[StellarYieldsMassiveStarData.Nm-1])){
    // use massive star data from PARSEC models (massive stars only! and wind only!)

    table = &StellarYieldsMassiveStarData;
  }
#endif

  return StellarYieldsGetYieldTablePosition(*table, i, j, M, metallicity);
}

int StellarYieldsGetPopIIIYieldTablePosition(int &i, const float &M){

    StellarYieldsDataType * table = &StellarYieldsPopIIIData;

    int j = -1;

    return StellarYieldsGetYieldTablePosition(*table, i, j, M, 0.0,
                                                    1 // extrapolate if mass is too big
                                                     );
}

float StellarYieldsInterpolatePopIIIYield(const int &i, const float &M, int atomic_number){

  StellarYieldsDataType * table = &StellarYieldsPopIIIData;

  float yield;

  float interp_M = M;

  if ((M < table->M[0]) || (M > table->M[table->Nm-1])){
    const float fudge = 0.000001; // mass needs to be within bounds
    interp_M  = (M < table->M[0]) ? table->M[0]*(1.0+fudge): table->M[table->Nm - 1] * (1.0-fudge);
  }


  float t,u;
  t = (interp_M - table->M[i]) / (table->M[i+1] - table->M[i]);

  float ll, ul; // lower left, lower right, upper right, upper left points
  if (atomic_number > 0) { // interpolate yields given atomic number

    int yield_index;
    yield_index     = GetYieldIndex(table->Ny, atomic_number);

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,0,yield_index,table->Nm,table->Nz);
    ll = table->Yields[index];
    ul = table->Yields[index + table->dm];
#else
    ll = table->Yields[i  ][0][yield_index];
    ul = table->Yields[i+1][0][yield_index];
#endif

  } else if (atomic_number == 0){ // interpolate total metal mass

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,0,0,table->Nm,table->Nz);
    ll = table->Metal_Mtot[index];
    ul = table->Metal_Mtot[index + table->dm];
#else
    ll = table->Metal_Mtot[i  ][0];
    ul = table->Metal_Mtot[i+1][0];
#endif

  } else if (atomic_number < 0 ){ // interpolate total mass

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,0,0,table->Nm,table->Nz);
    ll = table->Mtot[index];
    ul = table->Mtot[index + table->dm];
#else
    ll = table->Mtot[i  ][0];
    ul = table->Mtot[i+1][0];
#endif

  }

  return ((1.0 - t) * (ll)   +
          (      t) * (ul) ) * M / interp_M ;
}

float StellarYieldsInterpolateYield(int yield_type,
                                    const int &i, const int &j,
                                    const float &M, const float &metallicity, int atomic_number){
/* ------------------------------------------------------------------------
 * StellarYieldsInterpolateYield
 * ------------------------------------------------------------------------
 * Interpolation function for when the position in the table is already known
 * (i and j). Does some out of bounds checking just in case.
 *
 * yield_type switches between SN and wind yield tables, atomic_number
 * is the atomic number of the desired yield. If atomic_number < 0, the
 * total yield mass is returned. If atomic_number == 0, the total metal mass
 * is returned.
 * -------------------------------------------------------------------------*/

  StellarYieldsDataType * table;

  if (yield_type == 0){            // do supernova / end of life yields

    table = &StellarYieldsSNData;

    if( M > table->M[table->Nm - 1]){

      ENZO_FAIL("StellarYieldsInterpolateYield: No yields available for massive stars. Assuming all do direct collapse.");
    }
  } else if (yield_type == 1){     // do stellar wind yields

#ifdef NEWYIELDTABLES
    if (StellarYieldsAGBData.M == NULL){
      // use the old method
      if (M > StellarYieldsWindData.M[StellarYieldsWindData.Nm-1]){
        table = &StellarYieldsMassiveStarData;
      } else {
        table = &StellarYieldsWindData;
      }
    } else {
      /* New method picking between AGB and massive star */

      if (M > IndividualStarAGBThreshold){
        table = &StellarYieldsWindData;
      } else {
        table = &StellarYieldsAGBData;
      }

    } // end AGB == NULL

#else
    /* use NuGrid wind data if star is on grid, else use PARSEC massive star winds */

    if ( M > StellarYieldsWindData.M[StellarYieldsWindData.Nm - 1]){
      table = &StellarYieldsMassiveStarData;
    } else{
      table = &StellarYieldsWindData;
    }
#endif

  } else if (yield_type == 2){
    table = &StellarYieldsPopIIIData;
  }

  /* interpolate table */

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table->M[0]) || (M > table->M[table->Nm - 1])){
    if( IndividualStarExtrapolateYields ){
      const float fudge = 0.000001; // mass needs to be within bounds
      interp_M  = (M < table->M[0]) ? table->M[0]*(1.0+fudge): table->M[table->Nm - 1] * (1.0-fudge);

    } else{
       ENZO_FAIL("Interpolation mass off of yields table\n");
    }
  }

  if( (Z < table->Z[0]) || (Z > table->Z[table->Nz - 1])){
    if ( Z < table->Z[0] ){ // WARNING: see statement at top of file
      Z = table->Z[0];
    } else{
      ENZO_FAIL("Metallicity off of yields table\n");
    }
  }

  float t,u;
  t = (interp_M - table->M[i]) / (table->M[i+1] - table->M[i]);
  u = (Z  - table->Z[j]) / (table->Z[j+1] - table->Z[j]);

  if( (t<0) || (u<0)  || (t>1) || (u>1)){
      printf("Stellar yields interpolation issue - t = %"FSYM" u = %"FSYM"\n",t,u);
      printf("i, j, interp_M, Z %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n", i, j, interp_M, Z, table->M[i], table->Z[j]);
  }

  float ll, lr, ur, ul; // lower left, lower right, upper right, upper left points
  if (atomic_number > 0) { // interpolate yields given atomic number

    int yield_index;
    yield_index     = GetYieldIndex(table->Ny, atomic_number);

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,j,yield_index,table->Nm,table->Nz);

    ll = table->Yields[index];
    lr = table->Yields[index + table->dz];
    ur = table->Yields[index + table->dm + table->dz];
    ul = table->Yields[index + table->dm];

#else
    ll = table->Yields[i  ][j  ][yield_index];
    lr = table->Yields[i  ][j+1][yield_index];
    ur = table->Yields[i+1][j+1][yield_index];
    ul = table->Yields[i+1][j  ][yield_index];
#endif

  } else if (atomic_number == 0){ // interpolate total metal mass

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,j,0,table->Nm,table->Nz);

    ll = table->Metal_Mtot[index];
    lr = table->Metal_Mtot[index + table->dz];
    ur = table->Metal_Mtot[index + table->dm + table->dz];
    ul = table->Metal_Mtot[index + table->dm];

#else
    ll = table->Metal_Mtot[i  ][j  ];
    lr = table->Metal_Mtot[i  ][j+1];
    ur = table->Metal_Mtot[i+1][j+1];
    ul = table->Metal_Mtot[i+1][j  ];
#endif

  } else if (atomic_number < 0 ){ // interpolate total mass

#ifdef NEWYIELDTABLES
    int index = YIELD_INDEX(i,j,0,table->Nm,table->Nz);

    ll = table->Mtot[index];
    lr = table->Mtot[index + table->dz];
    ur = table->Mtot[index + table->dm + table->dz];
    ul = table->Mtot[index + table->dm];

#else
    ll = table->Mtot[i  ][j  ];
    lr = table->Mtot[i  ][j+1];
    ur = table->Mtot[i+1][j+1];
    ul = table->Mtot[i+1][j  ];
#endif

  }

  return ((1.0 - t)*(1.0 - u) * (ll)   +
          (1.0 - t)*(      u) * (lr)   +
          (      t)*(      u) * (ur)   +
          (      t)*(1.0 - u) * (ul) ) * M / interp_M ;
}

float StellarYieldsInterpolateYield(int yield_type,
                                    const float &M, const float &metallicity, int atomic_number){
  /* --------------------------------------------------------------------------
   * StellarYieldsInterpolateYield
   * --------------------------------------------------------------------------
   * Interpolate stellar yields table over mass and metallicity for given element
   * atomic number. If atomic number is negative, then instead interpolate along
   * the total ejected mass
   * ---------------------------------------------------------------------------- */

  StellarYieldsDataType * table;

#ifdef NEWYIELDTABLES
  if (M <= IndividualStarAGBThreshold){

    if (StellarYieldsAGBData.M == NULL){
      //  This is done for backwards compatability
      //  In original tables, winds contained AGB stars
      table = &StellarYieldsWindData;

    } else {
      // otherwise we are doig the new style
      table = &StellarYieldsAGBData;
    }

  } else{

    if (StellarYieldsMassiveStarData.M != NULL){
      // this is the OLD style for backwards compatability
      if (M <= StellarYieldsSNData.M[StellarYieldsSNData.Nm-1]){
        table = &StellarYieldsSNData;
      } else if (( M >= StellarYieldsMassiveStarData.M[0]) &&
                 ( M <= StellarYieldsMassiveStarData.M[StellarYieldsMassiveStarData.Nm-1])){
         // use massive star data from PARSEC models (massive stars only! and wind only!)

         table = &StellarYieldsMassiveStarData;
       }

    } else {
      // in the new model, both the SN table and the winds tables have
      // the same mass bins. In this case, it doesn't matter which we pick
      table = &StellarYieldsWindData;
    }
  }

# else
  if (yield_type == 0){            // do supernova / end of life yields

    table = &StellarYieldsSNData;

    if( M > table->M[table->Nm - 1]){
      ENZO_FAIL("StellarYieldsInterpolateYield: No yields available for massive stars. Assuming all do direct collapse.");
    }
  } else if (yield_type == 1){     // do stellar wind yields

    /* use NuGrid wind data if star is on grid, else use PARSEC massive star winds */
    if ( M > StellarYieldsWindData.M[StellarYieldsWindData.Nm - 1]){
      table = &StellarYieldsMassiveStarData;
    } else{
      table = &StellarYieldsWindData;
    }

  } // another if for SN1a
#endif
  /* interpolate table */

  int i, j;
  int width, bin_number;

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table->M[0]) || (M > table->M[table->Nm - 1])){
    if( IndividualStarExtrapolateYields ){
      // scale to most massive star
      const float fudge = 0.000001; // mass needs to be within bounds
      interp_M  = (M < table->M[0]) ? table->M[0]*(1.0+fudge): table->M[table->Nm - 1] * (1.0-fudge);
    } else{
      printf("StellarYieldsInterpolateYield: Mass out of bounds\n");
      printf("M = %"ESYM" for minimum M = %"ESYM" and maximum M = %"ESYM"\n", M, table->M[0], table->M[table->Nm-1]);
      return FAIL;
    }
  }

  if( (Z < table->Z[0]) || (Z > table->Z[table->Nz - 1])){

    if ( Z < table->Z[0] ){ // WARNING: see statement at top of file
      Z = table->Z[0];
    } else {

      printf("StellarYieldsInterpolateYield: Metallicity out of bounds\n");
      printf("Z = %"ESYM" for minimum Z = %"ESYM" and maximum Z = %"ESYM"\n", Z, table->Z[0], table->Z[table->Nz-1]);

      return FAIL;
    }
  }

  i = search_lower_bound(table->M, interp_M, 0, table->Nm, table->Nm);
  j = search_lower_bound(table->Z, Z, 0, table->Nz, table->Nz);

  return StellarYieldsInterpolateYield(yield_type,i,j,M,Z,atomic_number);

}

int GetYieldIndex(const int &number_of_yields, const int &Z){
/* --------------------------------------------------------
 * GetYieldIndex
 * --------------------------------------------------------
 * Unfortunate function, but way to map atomic number of
 * desired element species to the index where it is stored
 * in the yield data array. At first attemped to use
 * std::map for this purpose, but was having trouble getting
 * std::map to work in StellarYieldData.h
 * -------------------------------------------------------- */

 int i = 0;
 int index = -1;

 while( i < number_of_yields){
   if( StellarYieldsAtomicNumbers[i] == Z){
     index = i;
     break;
   }
   i++;
 }

  return index;
}
