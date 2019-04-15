/***********************************************************************
/
/  DETERMINE IF TWO BOXES OVERLAP
/
/  written by: Stephen Skory
/  date:       September, 2012
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"

int one_overlap(FLOAT a[MAX_DIMENSION], FLOAT b[MAX_DIMENSION],
                FLOAT s[MAX_DIMENSION], FLOAT t[MAX_DIMENSION],
                int dims) {
  int dim;
  int overlap = SUCCESS;
  for (dim = 0; dim < dims; dim++) {
    if ((a[dim] >= t[dim]) || (b[dim] <= s[dim])) {
      overlap = FAIL;
      break;
    }
  }
  if (overlap) {return SUCCESS;}
  else {return FAIL;}
}

int check_overlap(FLOAT a[MAX_DIMENSION], FLOAT b[MAX_DIMENSION],
                  FLOAT s[MAX_DIMENSION], FLOAT t[MAX_DIMENSION],
                  int dims, FLOAT period[MAX_DIMENSION],
                  int *shift, int skipto) {
  // a (left), b (right) - corners of box 1
  // s, t - corners of box 2
  // skipto - where to start the loops below.
  FLOAT a_temp[MAX_DIMENSION], b_temp[MAX_DIMENSION];
  int this_shift[MAX_DIMENSION];
  int shift0, shift1, shift2, dim, overlap, max1, max2, i1, i2, count = -1;
  // Test the simplest case, where they already overlap.
  //overlap = one_overlap(a, b, s, t, dims);
  //if (overlap) return SUCCESS;
  // If we're here, we need to test all cases.
  // Keep checking until we've exhausted all cases (as shift is handed back
  // to here from above).
  // This is to ensure that the box is shifted in all the dimensions it needs
  // to for reproducibility, meaning if we shift box 1 by an amount, we would
  // (swapping box 1 and box 2 in the algorithm) shift box 2 by the exact
  // opposite amount.
  // We shift box 1 around, keeping grid 2 static.
  // Decide how many dimensions we move in.
  max1 = (dims > 1) ? 2 : -1;
  max2 = (dims > 2) ? 2 : -1;
  i1 = (dims > 1) ? 1 : 0;
  i2 = (dims > 2) ? 2 : 0;
  // we start each loop at the *current* setting of shift.
  for (shift0 = -1; shift0 < 2; shift0++) {
    this_shift[0] = shift0;
    for (shift1 = -1; shift1 < max1; shift1++) {
      if (i2) this_shift[1] = shift1;
      for (shift2 = -1; shift2 < max2; shift2++) {
        if (i1) this_shift[2] = shift2;
        count += 1;
        if (count < skipto) continue;
        // if we're here, this is a new inspection.
        for (dim = 0; dim < dims; dim++) {
          a_temp[dim] = a[dim] + this_shift[dim] * period[dim];
          b_temp[dim] = b[dim] + this_shift[dim] * period[dim];
        }
        overlap = one_overlap(a_temp, b_temp, s, t, dims);
        if (overlap) {
          for (dim = 0; dim < dims; dim++) {
            shift[dim] = this_shift[dim];
            }
          return (count + 1);
         } // if overlap
      } // shift2
    } // shift1
  } // shift0
  
  // If we get here, they don't overlap.
  return -1;
}

