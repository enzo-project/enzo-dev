/*************************************************************************
/
/ STELLAR YIELDS DATA: CHEMICAL YIELDS, MECHANICAL LUMINOSITY, MASS EJECTA
/
/ Written by: A. Emerick
/ Date      : 3/6/16
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

float StellarYieldsInterpolateYield(StellarYieldsDataType table,
                                  float *M, float *metallicity, int atomic_number){

  /* interpolate table */

  int i, j;
  int width, bin_number;

  float Z = (*metallicity); // not in solar units

  if( (*M < table.M[0]) || (*M > table.M[table.NumberOfMassBins - 1])){
    printf("StellarYieldsInterpolateYield: Mass out of bounds\n");
    printf("M = %"ESYM" for minimum M = %"ESYM" and maximum M = %"ESYM"\n", *M, table.M[0], table.M[table.NumberOfMassBins-1]);
    return FAIL;
  }

  if( (Z < table.Z[0]) || (Z > table.Z[table.NumberOfMetallicityBins - 1])){
    printf("StellarYieldsInterpolateYield: Metallicity out of bounds\n");
    printf("Z = %"ESYM" for minimum Z = %"ESYM" and maximum Z = %"ESYM"\n", Z, table.Z[0], table.Z[table.NumberOfMetallicityBins-1]);
    return FAIL;
  }

  /* binary search over star mass */
  width = table.NumberOfMassBins / 2;
  i     = table.NumberOfMassBins / 2;

  while (width > 1){

    width /= 2;
    if ( *M > table.M[i])
      i += width;
    else if (*M < table.M[i])
      i -= width;
    else
      break;
  } // mass binary search

  if ( *M < table.M[i] ) i--;
  if ( *M < table.M[i] ) i--;

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

  t = (*M - table.M[i]) / (table.M[i+1] - table.M[i]);
  u = (Z  - table.Z[j]) / (table.Z[j+1] - table.Z[j]);

  int yield_index;
  yield_index     = GetYieldIndex(table.atomic_number, table.NumberOfYields, atomic_number);

  return (1.0 - t)*(1.0 - u) * (table.Yields[i][j][yield_index]) +
         (1.0 - t)*(      u) * (table.Yields[i][j][yield_index]) +
         (      t)*(      u) * (table.Yields[i][j][yield_index]) +
         (      t)*(1.0 - u) * (table.Yields[i][j][yield_index]);
}





int GetYieldIndex(int *atomic_numbers, int number_of_yields, int Z){
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
   if( atomic_numbers[i] == Z){
     index = i;
     break;
   }
 }

  return index;
}
