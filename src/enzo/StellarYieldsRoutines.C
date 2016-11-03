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
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "StellarYieldsRoutines.h"
#include "StarParticleData.h"

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

int StellarYieldsGetYieldTablePosition(int &i, int &j,
                                        const float &M, const float &metallicity){
  /* interpolate table */
  StellarYieldsDataType table = StellarYieldsSNData; // either is fine

  int width, bin_number;

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table.M[0]) || (M > table.M[table.NumberOfMassBins - 1])){
    if( IndividualStarExtrapolateYields ){
      // scale to most massive star
      interp_M  = table.M[table.NumberOfMassBins - 1] * 0.9999999;
    } else{
      printf("StellarYieldsInterpolateYield: Mass out of bounds\n");
      printf("M = %"ESYM" for minimum M = %"ESYM" and maximum M = %"ESYM"\n", M, table.M[0], table.M[table.NumberOfMassBins-1]);
      return FAIL;
    }
  }

  if( (Z < table.Z[0]) || (Z > table.Z[table.NumberOfMetallicityBins - 1])){

    if ( Z < table.Z[0] ){ // WARNING: see statement at top of file
      Z = table.Z[0];
    } else {

      printf("StellarYieldsInterpolateYield: Metallicity out of bounds\n");
      printf("Z = %"ESYM" for minimum Z = %"ESYM" and maximum Z = %"ESYM"\n", Z, table.Z[0], table.Z[table.NumberOfMetallicityBins-1]);

      return FAIL;
    }
  }

  /* binary search over star mass */
  width = table.NumberOfMassBins / 2;
  i     = table.NumberOfMassBins / 2;

  while (width > 1){

    width /= 2;
    if ( interp_M > table.M[i])
      i += width;
    else if (interp_M < table.M[i])
      i -= width;
    else
      break;
  } // mass binary search

  if ( interp_M < table.M[i] ) i--;
  if ( interp_M < table.M[i] ) i--;

  width = table.NumberOfMetallicityBins / 2;
  j     = table.NumberOfMetallicityBins / 2;

  while (width > 1){
    width /= 2;
    if ( Z > table.Z[j])
      j += width;
    else if ( Z < table.Z[j])
      j -= width;
    else
      break;
  } // metallicity binary search

  // search finds nearest bin - interpolation requires floored nearest bin
  if ( Z < table.Z[j]) j--;
  if ( Z < table.Z[j]) j--;

  return SUCCESS;
}

float StellarYieldsInterpolateYield(int yield_type,
                                    const int &i, const int &j,
                                    const float &M, const float &metallicity, int atomic_number){

  StellarYieldsDataType table;

  if (yield_type == 0){            // do supernova / end of life yields

    table = StellarYieldsSNData;

  } else if (yield_type == 1){     // do stellar wind yields

    table = StellarYieldsWindData;

  } // another if for SN1a

  /* interpolate table */

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table.M[0]) || (M > table.M[table.NumberOfMassBins - 1])){
    if( IndividualStarExtrapolateYields ){
      // scale to most massive star
      interp_M  = table.M[table.NumberOfMassBins - 1] * 0.9999999;
    }
  }

  if( (Z < table.Z[0]) || (Z > table.Z[table.NumberOfMetallicityBins - 1])){
    if ( Z < table.Z[0] ){ // WARNING: see statement at top of file
      Z = table.Z[0];
    }
  }

  float t,u;
  t = (interp_M - table.M[i]) / (table.M[i+1] - table.M[i]);
  u = (Z  - table.Z[j]) / (table.Z[j+1] - table.Z[j]);

  if( (t<0) || (u<0) ){
      printf("Stellar yields interpolation issue - t = %"FSYM" u = %"FSYM"\n",t,u);
      printf("i, j, interp_M, Z %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n", i, j, interp_M, Z, table.M[i], table.Z[j]);
  }

  float ll, lr, ur, ul; // lower left, lower right, upper right, upper left points
  if (atomic_number > 0) { // interpolate yields given atomic number

    int yield_index;
    yield_index     = GetYieldIndex(table.NumberOfYields, atomic_number);

    ll = table.Yields[i  ][j  ][yield_index];
    lr = table.Yields[i  ][j+1][yield_index];
    ur = table.Yields[i+1][j+1][yield_index];
    ul = table.Yields[i+1][j  ][yield_index];

  } else if (atomic_number == 0){ // interpolate total metal mass
    ll = table.Metal_Mtot[i  ][j  ];
    lr = table.Metal_Mtot[i  ][j+1];
    ur = table.Metal_Mtot[i+1][j+1];
    ul = table.Metal_Mtot[i+1][j  ];

  } else if (atomic_number < 0 ){ // interpolate total mass

    ll = table.Mtot[i  ][j  ];
    lr = table.Mtot[i  ][j+1];
    ur = table.Mtot[i+1][j+1];
    ul = table.Mtot[i+1][j  ];
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

  StellarYieldsDataType table;

  if (yield_type == 0){            // do supernova / end of life yields

    table = StellarYieldsSNData;

  } else if (yield_type == 1){     // do stellar wind yields

    table = StellarYieldsWindData;

  } // another if for SN1a

  /* interpolate table */

  int i, j;
  int width, bin_number;

  float Z = (metallicity); // not in solar units

  float interp_M = M;

  if( (M < table.M[0]) || (M > table.M[table.NumberOfMassBins - 1])){
    if( IndividualStarExtrapolateYields ){
      // scale to most massive star
      interp_M  = table.M[table.NumberOfMassBins - 1] * 0.9999999;
    } else{
      printf("StellarYieldsInterpolateYield: Mass out of bounds\n");
      printf("M = %"ESYM" for minimum M = %"ESYM" and maximum M = %"ESYM"\n", M, table.M[0], table.M[table.NumberOfMassBins-1]);
      return FAIL;
    }
  }

  if( (Z < table.Z[0]) || (Z > table.Z[table.NumberOfMetallicityBins - 1])){

    if ( Z < table.Z[0] ){ // WARNING: see statement at top of file
      Z = table.Z[0];
    } else {

      printf("StellarYieldsInterpolateYield: Metallicity out of bounds\n");
      printf("Z = %"ESYM" for minimum Z = %"ESYM" and maximum Z = %"ESYM"\n", Z, table.Z[0], table.Z[table.NumberOfMetallicityBins-1]);

      return FAIL;
    }
  }

  /* binary search over star mass */
  width = table.NumberOfMassBins / 2;
  i     = table.NumberOfMassBins / 2;

  while (width > 1){

    width /= 2;
    if ( interp_M > table.M[i])
      i += width;
    else if (interp_M < table.M[i])
      i -= width;
    else
      break;
  } // mass binary search

  if ( interp_M < table.M[i] ) i--;
  if ( interp_M < table.M[i] ) i--;

  width = table.NumberOfMetallicityBins / 2;
  j     = table.NumberOfMetallicityBins / 2;

  while (width > 1){
    width /= 2;
    if ( Z > table.Z[j])
      j += width;
    else if ( Z < table.Z[j])
      j -= width;
    else
      break;
  } // metallicity binary search

  // search finds nearest bin - interpolation requires floored nearest bin
  if ( Z < table.Z[j]) j--;
  if ( Z < table.Z[j]) j--;

  float t, u;

  t = (interp_M - table.M[i]) / (table.M[i+1] - table.M[i]);
  u = (Z  - table.Z[j]) / (table.Z[j+1] - table.Z[j]);

  float ll, lr, ur, ul; // lower left, lower right, upper right, upper left points
  if (atomic_number > 0) { // interpolate yields given atomic number

    int yield_index;
    yield_index     = GetYieldIndex(table.NumberOfYields, atomic_number);

    ll = table.Yields[i  ][j  ][yield_index];
    lr = table.Yields[i  ][j+1][yield_index];
    ur = table.Yields[i+1][j+1][yield_index];
    ul = table.Yields[i+1][j  ][yield_index];

  } else if (atomic_number == 0){ // interpolate total metal mass
    ll = table.Metal_Mtot[i  ][j  ];
    lr = table.Metal_Mtot[i  ][j+1];
    ur = table.Metal_Mtot[i+1][j+1];
    ul = table.Metal_Mtot[i+1][j  ];

  } else if (atomic_number < 0 ){ // interpolate total mass

    ll = table.Mtot[i  ][j  ];
    lr = table.Mtot[i  ][j+1];
    ur = table.Mtot[i+1][j+1];
    ul = table.Mtot[i+1][j  ];
  }

//  printf("t u ll lr ur ul M im %"ESYM" %"ESYM" %"ESYM " %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",t,u,ll,lr,ur,ul,M,interp_M);
  return ((1.0 - t)*(1.0 - u) * (ll)   +
          (1.0 - t)*(      u) * (lr)   +
          (      t)*(      u) * (ur)   +
          (      t)*(1.0 - u) * (ul) ) * M / interp_M ;
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
