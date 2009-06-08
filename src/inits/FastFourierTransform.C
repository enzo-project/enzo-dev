/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM (FORWARD OR INVERSE)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:  Robert Harkness
/  date:       May, 2003
/
/  PURPOSE:
/
/  INPUTS:
/      buffer - field to be FFTed
/      Rank   - rank of FFT
/      DimensionReal[] - declared dimensions of buffer
/      Dimension[]     - active dimensions of buffer
/      direction       - +1 forward, -1 inverse
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
 
#ifdef GOT_FFT
# undef GOT_FFT
#endif
 
//  Function prototypes
 
int FastFourierTransformSGIMATH(FLOAT *buffer, int Rank, int DimensionReal[],
                                int Dimension[], int direction, int type);
int FastFourierTransformSUNMATH(FLOAT *buffer, int Rank, int DimensionReal[],
                                int Dimension[], int direction, int type);
int FastFourierTransformIBMMATH(FLOAT *buffer, int Rank, int DimensionReal[],
                                int Dimension[], int direction, int type);
 
int FastFourierTransformPrepareComplex(FLOAT *buffer, int Rank, int DimensionReal[],
                                     int Dimension[], int direction, int type);
 
 
 
 
int FastFourierTransform(FLOAT *buffer, int Rank, int DimensionReal[],
			 int Dimension[], int direction, int type)
{
 
#if defined(IRIS4) && defined(SGI_MATH)
 
  /* Use SGI's library routines. */
 
  if (FastFourierTransformSGIMATH(buffer, Rank, DimensionReal,
				  Dimension, direction, type) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformSGIMATH.\n");
    return FAIL;
  }
 
#define GOT_FFT
#endif /* IRIS4 && SGI_MATH */
 
 
#if defined(SUN_MATH)
 
  /* Use SUN's library routines. */
 
#define GOT_FFT
#endif /* SUN_MATH */
 
 
#if defined(IBM_MATH)
 
  /* Use IBM's library routines. */
 
#define GOT_FFT
#endif /* IBM_MATH */
 
 
#ifndef GOT_FFT
 
  /* Catchall: if there is no library FFT, use Fortran Complex by Singleton */
 
  if (FastFourierTransformPrepareComplex(buffer, Rank, DimensionReal,
				         Dimension, direction, type) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransformPrepareComplex.\n");
    return FAIL;
  }
 
#endif /* GOT_FFT */
 
  return SUCCESS;
}
