/***********************************************************************
/
/  COMPUTE A PARALLEL FFT
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
 
int CommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder);
int FastFourierTransform(float *buffer, int Rank, int DimensionReal[],
			 int Dimension[], int direction, int type);
void PrintMemoryUsage(char *str);
 
int CommunicationParallelFFT(region *InRegion, int NumberOfInRegions,
			     region **OutRegion, int *NumberOfOutRegions,
			     int DomainDim[], int Rank,
			     int direction, int TransposeOnCompletion)
{
 
  /* Definitions. */
 
  int i, j, k, dim, size;
  float x, DomainCellSize[MAX_DIMENSION];

  PrintMemoryUsage("Enter FFT");

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    DomainCellSize[dim] = (DomainRightEdge[dim]-DomainLeftEdge[dim])/float(DomainDim[dim]);
 
  /* -------------------------------------- */
  /* Generate a description striped along dim 0. */
 
  region *strip0 = new region[NumberOfProcessors];
  for (i = 0, x = 0; i < NumberOfProcessors; i++) {
    strip0[i].StartIndex[0] = nint(x/(2.0*DomainCellSize[0]))*2;
    x += (DomainRightEdge[0]-DomainLeftEdge[0])/float(NumberOfProcessors);
    strip0[i].RegionDim[0] = nint(x/(2.0*DomainCellSize[0]))*2 -
      strip0[i].StartIndex[0];
    for (dim = 1; dim < MAX_DIMENSION; dim++) {
      strip0[i].StartIndex[dim] = 0;
      strip0[i].RegionDim[dim] = DomainDim[dim];
    }
    strip0[i].Processor = i;
    strip0[i].Data = NULL;
  }
 
  /* -------------------------------------- */
  /* Create regions striped along last dim. */
 
  region *strip1 = new region[NumberOfProcessors];
  int LastIndex = max(Rank-1, 1);
  for (i = 0, x = 0; i < NumberOfProcessors; i++) {
    strip1[i].StartIndex[LastIndex] = nint(x/DomainCellSize[LastIndex]);
    x += (DomainRightEdge[LastIndex] - DomainLeftEdge[LastIndex])/
      float(NumberOfProcessors);
    strip1[i].RegionDim[LastIndex] = nint(x/DomainCellSize[LastIndex]) -
                                   strip1[i].StartIndex[LastIndex];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (dim != LastIndex) {
	strip1[i].StartIndex[dim] = 0;
	strip1[i].RegionDim[dim] = DomainDim[dim];
      }
    strip1[i].Processor = i;
    strip1[i].Data = NULL;
  }
 
  /* -------------------------------------- */
 
  if (direction == FFT_FORWARD) {
 
    /* Copy initial regions to striped1 regions. */
 
//    fprintf(stderr, "FFT(%"ISYM"): initial -> strip1\n", MyProcessorNumber);
    if (CommunicationTranspose(InRegion, NumberOfInRegions, strip1,
			       NumberOfProcessors, NORMAL_ORDER) == FAIL) {
      ENZO_FAIL("Error in CommunicationTranspose.\n");
    }
 
    /* Compute the size of the each FFT 'slice', and the actual dimensions,
       here we assume that it is 2 less than the declared dimension. */
 
    int TempInts[] = {1,1,1};
    for (k = 0, size = 1; k < Rank; k++) {
      TempInts[k] = strip1[MyProcessorNumber].RegionDim[k] - ((k == 0)? 2:0);
      if (k != LastIndex)
	size *= strip1[MyProcessorNumber].RegionDim[k];
    }
 
    /* FFT all dims except last (real to complex). */
 
//    fprintf(stderr, "FFT(%"ISYM"): FFT strip1\n", MyProcessorNumber);
    for (i=0, j=0; i < strip1[MyProcessorNumber].RegionDim[LastIndex]; i++)
      if (strip1[MyProcessorNumber].Data != NULL) {
	if (FastFourierTransform(strip1[MyProcessorNumber].Data+j, LastIndex,
				 strip1[MyProcessorNumber].RegionDim, TempInts,
				 direction, REAL_TO_COMPLEX) == FAIL) {
	  ENZO_FAIL("Error in forward ParallelFFT call.\n");
	}
	j += size;
      } // end: loop over last dims
 
    if (Rank > 1) {
 
      /* Transpose to striped0 regions (reverse order within blocks). */
 
//      fprintf(stderr, "FFT(%"ISYM"): strip1 -> strip0\n", MyProcessorNumber);
      if (CommunicationTranspose(strip1, NumberOfProcessors, strip0,
			      NumberOfProcessors, TRANSPOSE_FORWARD) == FAIL) {
	ENZO_FAIL("Error in CommunicationTranspose (forward).\n");
      }
 
      /* FFT last dim (complex to complex).  Note that the data has been
         transposed so that the last dimension has the most rapidly changing
	 index. */
 
      int nffts = 1, fft_size = strip0[MyProcessorNumber].RegionDim[LastIndex];
      for (j = 0; j < Rank-1; j++)
	nffts *= strip0[MyProcessorNumber].RegionDim[j];
      nffts /= 2;  // since these are complex ffts
//      fprintf(stderr, "FFT(%"ISYM"): FFT strip0\n", MyProcessorNumber);
      for (j = 0; j < nffts; j++)
	if (FastFourierTransform(strip0[MyProcessorNumber].Data+j*fft_size*2,
				 1, &fft_size, &fft_size, direction,
				 COMPLEX_TO_COMPLEX) == FAIL) {
	  ENZO_FAIL("Error in forward ParallelFFT call.\n");
	}
 
    } // end: if (Rank > 1)
 
    /* Set the output to strip0 unless it's only 1D, then use strip1. */
 
    *OutRegion = ((Rank>1) ? strip0 : strip1);
    *NumberOfOutRegions = NumberOfProcessors;
 
    /* Return data to original layout if requested. */
 
    if (TransposeOnCompletion) {
//      fprintf(stderr, "FFT(%"ISYM"): Transpose on Completion\n", MyProcessorNumber);
      if (CommunicationTranspose(*OutRegion, NumberOfProcessors,
				InRegion, NumberOfInRegions,
				((Rank>1) ? TRANSPOSE_REVERSE : NORMAL_ORDER))
	  == FAIL) {
	ENZO_FAIL("Error in CommunicationTranspose.\n");
      }
      *OutRegion = InRegion;
      *NumberOfOutRegions = NumberOfInRegions;
    }
 
  } // end: if (direction == FFT_FORWARD)
 
  /* -------------------------------------- */
 
  if (direction == FFT_INVERSE) {
 
    /* Tranpose from original layout to strip0/1 (if necessary). */
 
    if (TransposeOnCompletion) {
//      fprintf(stderr, "FFT(%"ISYM"): Transpose on Start\n", MyProcessorNumber);
      if (CommunicationTranspose(InRegion, NumberOfInRegions,
			     ((Rank>1) ? strip0 : strip1), NumberOfProcessors,
			     ((Rank>1) ? TRANSPOSE_FORWARD : NORMAL_ORDER) )
	  == FAIL) {
	ENZO_FAIL("Error in CommunicationTranspose.\n");
      }
    } else {
 
      /* Otherwise the data has been preserved in OutRegion, use that. */
 
      if (Rank > 1) {
	delete [] strip0;
	strip0 = *OutRegion;
      } else {
	delete [] strip1;
	strip1 = *OutRegion;
      }
    }
 
    if (Rank > 1) {
 
      /* FFT last dim (complex to complex). */
 
      int nffts = 1, fft_size = strip0[MyProcessorNumber].RegionDim[LastIndex];
      for (j = 0; j < Rank-1; j++)
	nffts *= strip0[MyProcessorNumber].RegionDim[j];
      nffts /= 2; //  since these are complex ffts
      for (j = 0; j < nffts; j++)
	if (FastFourierTransform(strip0[MyProcessorNumber].Data+j*fft_size*2,
				 1, &fft_size, &fft_size, direction,
				 COMPLEX_TO_COMPLEX) == FAIL) {
	  ENZO_FAIL("Error in forward ParallelFFT call.\n");
	}
 
      /* Transpose to striped0 regions (reverse order within blocks). */
 
      if (CommunicationTranspose(strip0, NumberOfProcessors, strip1,
			      NumberOfProcessors, TRANSPOSE_REVERSE) == FAIL) {
	ENZO_FAIL("Error in CommunicationTranspose (reverse).\n");
      }
 
    } // end: if (Rank > 1)
 
    /* FFT all dims except last (real to complex). */
 
    int TempInts[] = {1,1,1};
    for (k = 0, size = 1; k < Rank; k++) {
      TempInts[k] = strip1[MyProcessorNumber].RegionDim[k] - ((k == 0)? 2:0);
      if (k != LastIndex)
	size *= strip1[MyProcessorNumber].RegionDim[k];
    }
    for (i=0, j=0; i < strip1[MyProcessorNumber].RegionDim[LastIndex]; i++)
      if (strip1[MyProcessorNumber].Data != NULL) {
	if (FastFourierTransform(strip1[MyProcessorNumber].Data+j, LastIndex,
				 strip1[MyProcessorNumber].RegionDim, TempInts,
				 direction, REAL_TO_COMPLEX) == FAIL) {
	  ENZO_FAIL("Error in forward ParallelFFT call.\n");
	}
	j += size;
      } // end: loop over last dims
 
    /* Copy initial regions to striped1 regions. */
 
    if (CommunicationTranspose(strip1, NumberOfProcessors,
			       InRegion, NumberOfInRegions,
			       NORMAL_ORDER) == FAIL) {
      ENZO_FAIL("Error in CommunicationTranspose.\n");
    }
    *OutRegion = InRegion;
    *NumberOfOutRegions = NumberOfProcessors;
 
  } // End: if (direction == FFT_INVERSE)
 
  /* Clean up. */
 
  if (*OutRegion != strip1)
    delete [] strip1;
  if (*OutRegion != strip0)

    delete [] strip0;

  PrintMemoryUsage("Exit FFT");
 
  return SUCCESS;
}
