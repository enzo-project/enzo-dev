/***********************************************************************
/
/  C WRAPPER FOR FORTRAN COMPLEX FAST FOURIER TRANSFORM
/
/  written by: Robert Harkness
/  date:       May, 2003
/  modified1:
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
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
// Function prototypes
 
extern "C" void FORTRAN_NAME(prefort)(float *data, float *work,
				      int *dim1, int *dim2, int *dim3,
				      int *rdim1, int *rdim2, int *rdim3,
				      int *rank, int *isign);
 
extern "C" void FORTRAN_NAME(fortfft)(float *data,
                                      int *dim1, int *dim2, int *dim3,
                                      int *ndim, int *isign);
 
 
 
 
int FastFourierTransformPrepareComplex(float *buffer, int Rank,
                                       int DimensionReal[],
                                       int Dimension[], int direction, int type)
{
 
  int dim, Dim[MAX_DIMENSION], DimReal[MAX_DIMENSION];
 
  // Error check.
 
  if (Rank < 1 || Rank > 3) {
    ENZO_VFAIL("Does not support Rank = %"ISYM"\n", Rank)
  }
 
  // Copy passed dims to Real dims to make sure they are at least 3d.
 
  for (dim = 0; dim < Rank; dim++) {
    Dim[dim] = Dimension[dim];
    DimReal[dim] = DimensionReal[dim];
  }
  for (dim = Rank; dim < MAX_DIMENSION; dim++) {
    Dim[dim] = 1;
    DimReal[dim] = 1;
  }
 
  // Compute size of required temporary buffer.
 
  int size = 1;
 
  for (dim = 0; dim < Rank; dim++)
    size *= Dimension[dim];
 
  /* a) real-to-complex or reverse.  This uses a complex-to-complex
        FFT so we must convert real part to complex; done in prefort */
 
  if (type == REAL_TO_COMPLEX) {
 
    // Allocate temporary buffer.
 
    float *tempbuffer = new float[2*size];
 
    // Prepare and carry out transform.
 
    direction *= -1;
 
    FORTRAN_NAME(prefort)(buffer, tempbuffer, &Dim[0], &Dim[1], &Dim[2],
                          &DimReal[0], &DimReal[1], &DimReal[2], &Rank,
                          &direction);
 
    direction *= -1;
 
  // Clean up.
 
    delete [] tempbuffer;
 
  }
 
  /* b) complex-to-complex.  Just pass the data to fortfft */
 
  if (type == COMPLEX_TO_COMPLEX) {
 
    direction *= -1;
 
    FORTRAN_NAME(fortfft)(buffer, &Rank, &Dim[0], &Dim[1], &Dim[2], &direction);
 
    direction *= -1;
 
   // Scale on inverse.
 
//    if (direction == FFT_INVERSE) {

//      float factor = 1.0/size;
//      for (int i = 0; i < size*2; i++)
//        buffer[i] *= factor;
//    }
 
  }
 
  return SUCCESS;
 
}
