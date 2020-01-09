/***********************************************************************
/
/  DETERMINE THE NEXT LARGEST AVAILABLE FFT SIZE
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE: This routine returns the next largest available size for
/           an FFT transform, given the library function available.
/
/         Supported platform     FFT library
/           SGI                    SGI_MATH
/           CONVEX                 VECLIB
/           ANY                    FOURN (from Numerical Recipes)
/
************************************************************************/

/* Some FFT packages support more than just radix-2 transforms.  This,
   however can use an excessive number of Green's functions.  The following
   define flag indicates that only radix-2 transforms should be used, even
   if other options are available. */

// #define USE_ONLY_RADIX2
#define LIMITED_RANGE

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"

#ifdef GOT_MACHINE
#undef GOT_MACHINE
#endif

int NextLargestFFTSize(int dimension)
{

  /* --------------------------------------------------------------------- */
  /* Define the available sizes for SGI's WITH SGI_MATH */

#if defined(IRIS4) && defined(SGI_MATH) && !defined(USE_ONLY_RADIX2)
 #define RADIX235
 #define GOT_MACHINE
#endif /* IRIS4 && SGI_MATH */

  /* --------------------------------------------------------------------- */
  /* Define the available sizes for CONVEX WITH VECLIB. */

#if (defined(CONVEX) || defined(SPP)) && defined(VECLIB) && !defined(USE_ONLY_RADIX2)
 #define RADIX235
 #define GOT_MACHINE
#endif /* (CONVEX || SPP) && VECLIB */

  /* --------------------------------------------------------------------- */
  /* Define the available sizes for CONVEX WITH VECLIB. */

#if defined(FFTW) && !defined(USE_ONLY_RADIX2)
 #define RADIX235
 #define GOT_MACHINE
#endif /* (CONVEX || SPP) && VECLIB */

  /* --------------------------------------------------------------------- */
  /* If the machine was not one of the above, use this. */

#ifndef GOT_MACHINE

  /* FOURN is power of 2. */

  static int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

#endif /* GOT_MACHINE */  

  /* If FFT radix 2,3,5.  Many values over 64 were skipped. */

#ifdef RADIX235
 #ifdef LIMITED_RANGE
  static int sizes[] = {2, 4, 8, 16, 24, 32, 40, 48, 64, 96, 128, 256, 512, 1024};
 #else /* LIMITED_RANGE */
  static int sizes[] = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20,
			24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60,
			64, 100, 128, 150, 200, 256, 400, 512, 1024};
 #endif /* LIMITED_RANGE */
#endif /* RADIX 235 */

  /* --------------------------------------------------------------------- */
  /* Error check. */

#define NUM_SIZES (sizeof sizes / sizeof sizes[0])

  if (dimension > sizes[NUM_SIZES-1]) {
    fprintf(stderr, "NextLargestFFTSize: %d too large!\n", dimension);
    exit(FAIL);
  }

  /* Loop through available sizes. */

  int i = 0;
  while (sizes[i] < dimension && i < NUM_SIZES-1)
    i++;

  return sizes[i];

}
