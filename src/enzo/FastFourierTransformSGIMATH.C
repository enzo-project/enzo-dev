/***********************************************************************
/
/  COMPUTE A FAST FOURIER TRANSFORM USING THE SGI_MATH LIBRARY
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/      buffer - field to be FFTed
/      Rank   - rank of FFT
/      DimensionReal[] - declared dimensions of buffer
/      Dimension[]     - active dimensions of buffer
/      dir             - +1 forward, -1 inverse
/      type            - RTOC - real-to-complex, CTOC - complex-to-complex
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
/* Defines */
 
#if defined(IRIS4) && defined(SGI_MATH)
#include <fft.h>
#define GOT_FFT
#endif /* IRIS4 && SGI_MATH */
 
/* Call correct routine for single/double precision. */
 
#ifdef CONFIG_BFLOAT_4
#define RTOC(X)          scfft##X
#define RTOC_INVERSE(X)  csfft##X
#define RTOC_SCALE(X)    sscal##X
#define CTOC(X)           cfft##X
#define CTOC_SCALE(X)    cscal##X
#define CMPLX_CAST    (complex *)
#endif /* r4 */
 
#ifdef CONFIG_BFLOAT_8
#define RTOC(X)          dzfft##X
#define RTOC_INVERSE(X)  zdfft##X
#define RTOC_SCALE(X)    dscal##X
#define CTOC(X)           zfft##X
#define CTOC_SCALE(X)    zscal##X
#define CMPLX_CAST    (zomplex *)
#endif /* r8 */
 
/* Start routine. */
 
int FastFourierTransformSGIMATH(float *buffer, int Rank, int DimensionReal[],
				int Dimension[], int dir, int type)
{
 
#if defined(IRIS4) && defined(SGI_MATH)
 
  /* Use SGI's routines.  Note: the direction flag is defined backwards to the
     usual sense (duh). */
 
  /* 1D transform. */
 
  if (Rank == 1) {
    float *Work = new float[4*(DimensionReal[0]+15)];
    if (type == REAL_TO_COMPLEX) {
      RTOC(1di)(Dimension[0], Work);
      if (dir == FFT_FORWARD)
	RTOC(1du)(dir*-1, Dimension[0], buffer, 1, Work);
      else {
	RTOC_INVERSE(1du)(dir*-1, Dimension[0], buffer, 1, Work);
	RTOC_SCALE(1d)(Dimension[0], 1.0/float(Dimension[0]), buffer, 1);
      }
    } else {
      CTOC(1di)(Dimension[0], CMPLX_CAST Work);
      CTOC(1d)(dir*-1, Dimension[0], CMPLX_CAST buffer, 1, CMPLX_CAST Work);
      if (dir == FFT_INVERSE)
	CTOC_SCALE(1d)(Dimension[0], 1.0/float(Dimension[0]),
		       CMPLX_CAST buffer, 1);
    }
    delete Work;
  } // end of 1D FFT
 
  /* 2D transform. */
 
  if (Rank == 2) {
    float *Work = new float[4*(DimensionReal[0]+15 + 2*DimensionReal[1]+15)];
    if (type == REAL_TO_COMPLEX) {
      RTOC(2dui)(Dimension[0], Dimension[1], Work);
      if (dir == FFT_FORWARD)
	RTOC(2du)(dir*-1, Dimension[0], Dimension[1], buffer,
		  DimensionReal[0], Work);
      else {
	RTOC_INVERSE(2du)(dir*-1, Dimension[0], Dimension[1], buffer,
		  DimensionReal[0], Work);
	RTOC_SCALE(2d)(Dimension[0], Dimension[1],
		       1.0/float(Dimension[0]*Dimension[1]), buffer,
		       DimensionReal[0]);
      }
    } else {
      CTOC(2di)(Dimension[0], Dimension[1], CMPLX_CAST Work);
      CTOC(2d)(dir*-1, Dimension[0], Dimension[1], CMPLX_CAST buffer,
		  DimensionReal[0], CMPLX_CAST Work);
      if (dir == FFT_INVERSE)
	CTOC_SCALE(2d)(Dimension[0], Dimension[1],
		       1.0/float(Dimension[0]*Dimension[1]), CMPLX_CAST buffer,
		       DimensionReal[0]);
    }
    delete Work;
  } // end of 2D FFT
 
  /* 3D transform. */
 
  if (Rank == 3) {
    float *Work = new float[4*(  DimensionReal[0]+15 + 2*DimensionReal[1]+15 +
			       2*DimensionReal[2]+15)];
    if (type == REAL_TO_COMPLEX) {
      RTOC(3dui)(Dimension[0], Dimension[1], Dimension[2], Work);
      if (dir == FFT_FORWARD)
	RTOC(3du)(dir*-1, Dimension[0], Dimension[1], Dimension[2],
		  buffer, DimensionReal[0], DimensionReal[1], Work);
      else {
	RTOC_INVERSE(3du)(dir*-1, Dimension[0], Dimension[1], Dimension[2],
		  buffer, DimensionReal[0], DimensionReal[1], Work);
	RTOC_SCALE(3d)(Dimension[0], Dimension[1], Dimension[2],
		       1.0/float(Dimension[0]*Dimension[1]*Dimension[2]),
		       buffer, DimensionReal[0], DimensionReal[1]);
      }
    } else {
      CTOC(3di)(Dimension[0], Dimension[1], Dimension[2], CMPLX_CAST Work);
      CTOC(3d)(dir*-1, Dimension[0], Dimension[1], Dimension[2],
	       CMPLX_CAST buffer, DimensionReal[0], DimensionReal[1],
	       CMPLX_CAST Work);
      if (dir == FFT_FORWARD)
	CTOC_SCALE(3d)(Dimension[0], Dimension[1], Dimension[2],
		       1.0/float(Dimension[0]*Dimension[1]*Dimension[2]),
		       CMPLX_CAST buffer, DimensionReal[0], DimensionReal[1]);
    }
    delete Work;
  } // end of 3D FFT
 
  return SUCCESS;
 
#else /* IRIS4 && SGI_MATH */
 
  /* This is an error. */
 
  ENZO_FAIL("What are we doing here!?!\n");

 
#endif /* IRIS4 && SGI_MATH */
 
}
