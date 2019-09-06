/***********************************************************************
/
/  GRID CLASS (COMPUTE THE ACCLERATION FIELDS FROM THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int NextLargestFFTSize(int dimension);
int FastFourierTransform(float *buffer, int Rank, int DimensionReal[],
			 int Dimension[], int direction, int type);
int ReportMemoryUsage(char *header = NULL);
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);

#define INDEX_ACCEL(a,b,c) ( ((c)*GridDimension[1] + (b))*GridDimension[0] + (a) )
#define INDEX_GRAV(a,b,c) ( ((c)*GravitatingMassFieldDimension[1] + (b))*GravitatingMassFieldDimension[0] + (a) )
#define INDEX_FFT(a,b,c) ( ((c)*FFTDimensionReal[1] + (b))*FFTDimensionReal[0] + (a) )


int grid::ComputePotentialFieldAPM(int RefinementFactor)
{
  /* declarations */

  int i, j, k, dim, Zero[MAX_DIMENSION], bufferindex, gravityindex;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Zero[dim] = 0;

  /* Compute adot/a at time = t+1/2dt (time-centered). */

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }

  /* Determine the smallest possible dimensions that are:
     (a) larger than or equal in size to the gravity grid and
     (b) can be directly transformed.  This last quantity is
         implementation dependant (see NextLargestFFTSize). */

  int FFTDimension[MAX_DIMENSION], FFTDimensionReal[MAX_DIMENSION];
  int FFTDimensionGrav[MAX_DIMENSION], FFTStartIndex[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    FFTDimensionGrav[dim] = GravitatingMassFieldDimension[dim];
    FFTDimension[dim] = GravitatingMassFieldDimension[dim];
    FFTStartIndex[dim] = 0;

    /* If this is a topgrid periodic, then set the parameters to copy out just
       the active region. */

    if (GravityBoundaryType == TopGridPeriodic) {
      FFTStartIndex[dim] = nint( (GridLeftEdge[dim] - GravitatingMassFieldLeftEdge[dim]) / GravitatingMassFieldCellSize);
      FFTDimension[dim] = nint( (GridRightEdge[dim]-GridLeftEdge[dim])/GravitatingMassFieldCellSize);
      FFTDimensionGrav[dim] = FFTDimension[dim];
      if (NumberOfProcessors > 1) {
	fprintf(stderr, "This gravity type doesn't support parallel yet.\n");
	return FAIL;
      }
    }

    /* If subgrid isolated, then add one convolution kernel room for zero-padding. */

    if (GravityBoundaryType == SubGridIsolated)
      FFTDimension[dim] += nint(RefinementFactor*S2ParticleSize) + FFT_SAFETY_FACTOR;

    /* Set FFT dimension to next largest size we can do. */

    FFTDimension[dim] = NextLargestFFTSize(FFTDimension[dim]);
    FFTDimensionReal[dim] = FFTDimension[dim];
  }

  /* Add two to the declared dimensions of the FFT buffer since some
     real-to-complex transforms require this. */

  FFTDimensionReal[0] += 2;

  /* Compute the size of the FFT buffer and the gravitating field.
     Also, If this is the top grid and it's periodic, check to make
     sure the FFT size is exactly the same size as GravitatingMassField. */

  int FFTSize = 1, GravSize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    FFTSize  *= FFTDimensionReal[dim];
    GravSize *= GravitatingMassFieldDimension[dim];
    if (GravityBoundaryType == TopGridPeriodic &&
	FFTDimension[dim] != FFTDimensionGrav[dim]) {
      fprintf(stderr, "Topgrid cannot be periodic with these dimensions.\n");
      fprintf(stderr, "  (see NextLargestFFTSize.cc)\n");
      return FAIL;
    }
  }

  /* Allocate and clear two temporary buffers for the FFT. */

  float *buffer1 = new float[FFTSize];
  float *buffer2 = new float[FFTSize];
  if (buffer1 == NULL || buffer2 == NULL) {
    fprintf(stderr, "Grid_ComputeAccelerationField: malloc error (out of memory?)\n");
    return FAIL;
  }
  for (i = 0; i < FFTSize; i++)
    buffer1[i] = 0.0;

  /* Copy the GravitatingMassField into this buffer (first, clear unused
     indexes) and multiply by the gravitational constant (and by the cell
     width squared because the Green's function is unitless). */

  float factor = GravitationalConstant*GravitatingMassFieldCellSize
                                      *GravitatingMassFieldCellSize/a;
  printf("factor = %g\n", factor);
  for (k = 0; k < FFTDimensionGrav[2]; k++) {
    for (j = 0; j < FFTDimensionGrav[1]; j++) {
      for (i = 0; i < FFTDimensionGrav[0]; i++) {
	buffer1[INDEX_FFT(i,j,k)] = GravitatingMassField[INDEX_GRAV(i+FFTStartIndex[0], j+FFTStartIndex[1], k+FFTStartIndex[2])] * factor;

      }
    }
  }

  /* Transform the buffer into the fourier domain (in place). */

  if (FastFourierTransform(buffer1, GridRank, FFTDimensionReal, FFTDimension,
			   FFT_FORWARD, REAL_TO_COMPLEX) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransform (forward).\n");
    return FAIL;
  }

  /* Compute appropriate Green's function (not including D(k)). */

  if (GreensFunction.PrepareGreensFunction(GridRank, FFTDimensionReal,
	       FFTDimension, GravityBoundaryType, RefinementFactor) == FAIL) {
    fprintf(stderr, "Error in GreensFunction->PrepareGreensFunction.\n");
    return FAIL;
  }

  /* Multiply the Green's function by the transformed gravity field
     and leave the result (the potential) in buffer1. */

  if (GreensFunction.MultiplyGreensFunction(GridRank, FFTDimensionReal,
		        FFTDimension, buffer1) == FAIL) {
    fprintf(stderr, "Error in GreensFunction->MultiplyGreensFunction.\n");
    return FAIL;
  }

  /* Compute the potential with inverse FFT */

  if (FastFourierTransform(buffer1, GridRank, FFTDimensionReal,
			   FFTDimension, FFT_INVERSE, REAL_TO_COMPLEX) == FAIL) {
    fprintf(stderr, "Error in FastFourierTransform (inverse).\n");
    return FAIL;
  }

  /* Allocate potential field and copy. */

  if (PotentialField == NULL)
    PotentialField = new float[GravSize];

  //printf("In Grid_ComputePotentialFiedlAPM.c, set phi to 7\n");

  float maxpot = 0.0;
  float minpot = 0.0;
  for (k = 0; k < FFTDimensionGrav[2]; k++) {
    for (j = 0; j < FFTDimensionGrav[1]; j++) {
      for (i = 0; i < FFTDimensionGrav[0]; i++) {
	PotentialField[INDEX_GRAV(i+FFTStartIndex[0], j+FFTStartIndex[1], k+FFTStartIndex[2])] = buffer1[INDEX_FFT(i,j,k)];
      }
    }
  }

  /* If requested, compute the gravitational potential sum now. */

  int PotStart[MAX_DIMENSION], PotEnd[MAX_DIMENSION];
  float CellVolume = 1;
  for (dim = 0; dim < GridRank; dim++) {
    PotStart[dim] = nint((GridLeftEdge[dim] -
			  GravitatingMassFieldLeftEdge[dim])/
			 GravitatingMassFieldCellSize);
    PotEnd[dim] = nint((GridRightEdge[dim] -
			GravitatingMassFieldLeftEdge[dim])/
		       GravitatingMassFieldCellSize) - 1;
    CellVolume *= GravitatingMassFieldCellSize;
  }

  printf("PotStart/End = %d %d %d / %d %d %d\n", PotStart[0], PotStart[1],
	 PotStart[2], PotEnd[0], PotEnd[1], PotEnd[2]);

  PotentialSum = 0;
  for (k = PotStart[2]; k <= PotEnd[2]; k++)
    for (j = PotStart[1]; j <= PotEnd[1]; j++) {
      gravityindex = k*FFTDimensionGrav[0]*FFTDimensionGrav[1] +
	j*FFTDimensionGrav[0] + PotStart[0];
      for (i = PotStart[0]; i <= PotEnd[0]; i++, gravityindex++)
	PotentialSum += 0.5*PotentialField[gravityindex] *
	  GravitatingMassField[gravityindex]*
	  CellVolume;
    }
  printf("PotentialSum = %g\n", PotentialSum);

  /* clean up */

  GreensFunction.Release();

  delete [] buffer1;
  delete [] buffer2;

  return SUCCESS;
}
