#ifndef UNIGRID_TRANSPOSE
/***********************************************************************
/
/  COPY A SET OF REGIONS BETWEEN PROCESSORS (TRANSPOSE)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
extern "C" void FORTRAN_NAME(copy3dft)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
extern "C" void FORTRAN_NAME(copy3drt)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
 
 
 
int CommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{
 
  /* Declarations. */
 
  int dim, n, i, j, size, index, Zero[] = {0,0,0};
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];
  float *ReceiveBuffer, *SendBuffer;
 
  //  fprintf(stderr, "CT(%"ISYM"): start From=%"ISYM"  To=%"ISYM"\n", MyProcessorNumber,
  //	  NumberOfFromRegions, NumberOfToRegions);
							
  region *Sends = new region[NumberOfFromRegions];
  region *Receives = new region[NumberOfToRegions];
 
  /* Loop over processor jumps (number of processors ahead to send). */
 
  for (n = 0; n < NumberOfProcessors; n++) {
 
    /* Copy regions into communication buffer (or just set buffer
       if there is only one FromRegion per processor). */
 
    int sends = 0, receives = 0;
    int SendSize = 0, ReceiveSize = 0;
 
    for (j = 0; j < NumberOfFromRegions; j++)
      for (i = 0; i < NumberOfToRegions; i++)
        if ((ToRegion[i].Processor - FromRegion[j].Processor +
             NumberOfProcessors) % NumberOfProcessors == n &&
	    (MyProcessorNumber == FromRegion[j].Processor ||
	     MyProcessorNumber ==   ToRegion[i].Processor)) {
 
	  /* Determine if there is an overlap. */
 
	  size = 1;
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    LeftIndex[dim] = max(FromRegion[j].StartIndex[dim],
                                   ToRegion[i].StartIndex[dim]);
	    RightIndex[dim] = min(
	     FromRegion[j].StartIndex[dim] + FromRegion[j].RegionDim[dim],
	       ToRegion[i].StartIndex[dim] +   ToRegion[i].RegionDim[dim])-1;
	    size *= max(RightIndex[dim] - LeftIndex[dim] + 1, 0);
	  }
 
	  /* If there is an overlap, add it to the list of sends/receives. */
 
	  if (size > 0) {
 
	    if (MyProcessorNumber == FromRegion[j].Processor) {
//	      fprintf(stderr, "CT(%"ISYM"): from: i,j=%"ISYM" %"ISYM"  %"ISYM"->%"ISYM"\n",
//		      MyProcessorNumber, j, i, FromRegion[j].Processor,
//		      ToRegion[i].Processor);
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		Sends[sends].StartIndex[dim] = LeftIndex[dim] -
		  FromRegion[j].StartIndex[dim];
		Sends[sends].RegionDim[dim] = RightIndex[dim] -
		  LeftIndex[dim] + 1;
	      }
	      SendSize += size;
	      Sends[sends++].Processor = j;
	    }
 
	    if (MyProcessorNumber == ToRegion[i].Processor) {
//	      fprintf(stderr, "CT(%"ISYM"): to: %"ISYM"->%"ISYM" from proc %"ISYM"\n",
//		      MyProcessorNumber, j, i, FromRegion[j].Processor);
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		Receives[receives].StartIndex[dim] = LeftIndex[dim] -
		  ToRegion[i].StartIndex[dim];
		Receives[receives].RegionDim[dim] = RightIndex[dim] -
		  LeftIndex[dim] + 1;
	      }
	      ReceiveSize += size;
	      Receives[receives++].Processor = i;
	    }
 
	  } // end: if (size > 0)
 
        } // end: if (proc jump == n)
 
    /* Allocate buffer and copy data into buffer. */
 
    ReceiveBuffer = NULL;
    SendBuffer = new float[SendSize];
 
    index = 0;
 
//  fprintf(stderr, "CT(%"ISYM"): sends = %"ISYM"  SendSize = %"ISYM"\n", MyProcessorNumber,
//	    sends, SendSize);
 
    for (i = 0; i < sends; i++) {
      j = Sends[i].Processor;
      if (FromRegion[j].Data == NULL) {
	int ijk, rsize = FromRegion[j].RegionDim[0]*FromRegion[j].RegionDim[1]*FromRegion[j].RegionDim[2];
	FromRegion[j].Data = new float[rsize]; // allocate FROM field if not already done
	for (ijk = 0; ijk < rsize; ijk++)
	  FromRegion[j].Data[ijk] = 0;
      }
      if (TransposeOrder == TRANSPOSE_REVERSE)
	FORTRAN_NAME(copy3drt)(FromRegion[j].Data, SendBuffer+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1,
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1,
			   Sends[i].StartIndex+2);
      else
	FORTRAN_NAME(copy3d)(FromRegion[j].Data, SendBuffer+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1,
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1,
			   Sends[i].StartIndex+2);
      index += Sends[i].RegionDim[0]*Sends[i].RegionDim[1]*
	       Sends[i].RegionDim[2];
    }
 
    /* shift buffer by n processors. */
 
    if (n > 0) {
 
      ReceiveBuffer = new float[ReceiveSize];
 
#ifdef USE_MPI
 
      MPI_Status status;
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg Count;
      MPI_Arg RecvCount;
      MPI_Arg Source;
      MPI_Arg Dest;
 
      int ToProc = (MyProcessorNumber + n) % NumberOfProcessors;
      int FromProc = (MyProcessorNumber - n + NumberOfProcessors) %
	NumberOfProcessors;
 
#ifdef MPI_INSTRUMENTATION
      starttime = MPI_Wtime();
#endif
 
//      fprintf(stderr, "CT(%"ISYM"): MPI SS/RS = %"ISYM"/%"ISYM" From/To = %"ISYM" %"ISYM"\n",
//	      MyProcessorNumber, SendSize, ReceiveSize, FromProc, ToProc);

      Count = SendSize;
      RecvCount = ReceiveSize;
      Source = FromProc;
      Dest = ToProc;
       
      if (MPI_Sendrecv((void*) SendBuffer, Count, DataType, Dest,
	       MPI_TRANSPOSE_TAG, (void*) ReceiveBuffer, RecvCount,
	       DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
	       &status) != MPI_SUCCESS) {
	fprintf(stderr, "Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber,
		status.MPI_ERROR);
	return FAIL;
      }
 
#ifdef MPI_INSTRUMENTATION
      endtime = MPI_Wtime();
      timer[14] += endtime-starttime;
      counter[14] ++;
      CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
    } else {
      ReceiveBuffer = SendBuffer;
      SendBuffer = NULL;
    }
 
    /* Copy from communication buffer back to regions. */
 
    index = 0;
 
//    fprintf(stderr, "CT(%"ISYM"): receives = %"ISYM"\n", MyProcessorNumber, receives);
 
    for (i = 0; i < receives; i++) {

      j = Receives[i].Processor;

#ifdef ISOLATED_GRAVITY
      if (ToRegion[j].Data == NULL) {
	int ijk, rsize = ToRegion[j].RegionDim[0]*ToRegion[j].RegionDim[1]*ToRegion[j].RegionDim[2];
	ToRegion[j].Data = new float[rsize];
	for (ijk = 0; ijk < rsize; ijk++)
	  ToRegion[j].Data[ijk] = 0;
      }
#else /* ISOLATED_GRAVITY */
      if (ToRegion[j].Data == NULL)
	ToRegion[j].Data = new float[ToRegion[j].RegionDim[0]*
	      ToRegion[j].RegionDim[1]*ToRegion[j].RegionDim[2]];
#endif /* ISOLATED_GRAVITY */

      if (TransposeOrder == TRANSPOSE_FORWARD)
	FORTRAN_NAME(copy3dft)(ReceiveBuffer+index, ToRegion[j].Data,
		       Receives[i].RegionDim, Receives[i].RegionDim+1,
			   Receives[i].RegionDim+2,
		       ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
			   ToRegion[j].RegionDim+2,
		       Receives[i].StartIndex, Receives[i].StartIndex+1,
			   Receives[i].StartIndex+2,
		       Zero, Zero+1, Zero+2);
      else
	FORTRAN_NAME(copy3d)(ReceiveBuffer+index, ToRegion[j].Data,
		       Receives[i].RegionDim, Receives[i].RegionDim+1,
			   Receives[i].RegionDim+2,
		       ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
			   ToRegion[j].RegionDim+2,
		       Receives[i].StartIndex, Receives[i].StartIndex+1,
			   Receives[i].StartIndex+2,
		       Zero, Zero+1, Zero+2);
      index += Receives[i].RegionDim[0]*Receives[i].RegionDim[1]*
	       Receives[i].RegionDim[2];
    }
 
    /* Clean up. */
 
//    fprintf(stderr, "CT(%"ISYM"): end jump %"ISYM"\n", MyProcessorNumber, n);
 
    delete [] SendBuffer;
    delete [] ReceiveBuffer;
 
  } // end: loop over processors jumps
 
  /* Delete FromRegion data. */
 
//  fprintf(stderr, "CT(%"ISYM"): Deleting FromRegions\n", MyProcessorNumber);
 
  for (i = 0; i < NumberOfFromRegions; i++) {
    delete [] FromRegion[i].Data;
    FromRegion[i].Data = NULL;
  }
 
  /* Clean up. */
 
  delete [] Sends;
  delete [] Receives;
 
  return SUCCESS;
}
#endif
