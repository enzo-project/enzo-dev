#ifdef UNIGRID_TRANSPOSE
 
/***********************************************************************
/
/  COPY A SET OF REGIONS BETWEEN PROCESSORS (TRANSPOSE)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:  Giri Chukkapalli & Robert Harkness
/  date:       August, 2004
/              Static book-keeping for unigrid only!
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
 
//  Global definition for static book-keeping
 
commSndRcv *cSndRcv;
 
 
 
 
int CommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{
 
  int dim, n, i, j, size, index, Zero[] = {0,0,0};
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];
  float *ReceiveBuffer, *SendBuffer;
 
  int sends, receives;
  int SendSize, ReceiveSize;
 
  // Declarations for static unigrid book-keeping
 
  struct region *Sends, *Receives;
  int oSet;
 
  // Offset for each call - there are 6 axes
 
  oSet = NumberOfProcessors*((First_Pass+6)%6);
 
  //  fprintf(stderr, "CT(%"ISYM"): start From=%"ISYM"  To=%"ISYM"  offset=%"ISYM"\n", MyProcessorNumber,
  //          NumberOfFromRegions, NumberOfToRegions);
 
  // The first time this routine is called in each of 6 axes,
  // collect the static routing information for unigrid
 
  // Create a static table of routing information, once only
 
  if(First_Pass == 0)
    cSndRcv  = new commSndRcv[6*NumberOfProcessors];	
 
  if(First_Pass == 0)
    fprintf(stderr, "CAUTION: Unigrid FFT routing is effective!\n");
 
  if(First_Pass < 6)  {
 
  // Loop over processor jumps (number of processors ahead to send)
 
  for (n = 0; n < NumberOfProcessors; n++) {
 
    cSndRcv[oSet+n].Sends = new region[NumberOfFromRegions];
    cSndRcv[oSet+n].Receives = new region[NumberOfToRegions];
 
    /* Copy regions into communication buffer (or just set buffer
       if there is only one FromRegion per processor) */
 
    sends = 0, receives = 0;
    SendSize = 0, ReceiveSize = 0;
 
    for (j = 0; j < NumberOfFromRegions; j++)
      for (i = 0; i < NumberOfToRegions; i++)
        if ((ToRegion[i].Processor - FromRegion[j].Processor +
             NumberOfProcessors) % NumberOfProcessors == n &&
	    (MyProcessorNumber == FromRegion[j].Processor ||
	     MyProcessorNumber ==   ToRegion[i].Processor)) {
 
	  // Determine if there is an overlap
 
	  size = 1;
 
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    LeftIndex[dim] = max(FromRegion[j].StartIndex[dim],
                                   ToRegion[i].StartIndex[dim]);
	    RightIndex[dim] = min(
	     FromRegion[j].StartIndex[dim] + FromRegion[j].RegionDim[dim],
	       ToRegion[i].StartIndex[dim] +   ToRegion[i].RegionDim[dim])-1;
	    size *= max(RightIndex[dim] - LeftIndex[dim] + 1, 0);
	  }
 
	  // If there is an overlap, add it to the list of sends/receives
 
	  if (size > 0) {
 
	    if (MyProcessorNumber == FromRegion[j].Processor) {
 
              //  fprintf(stderr, "CT(%"ISYM"): from: i,j=%"ISYM" %"ISYM"  %"ISYM"->%"ISYM"\n",
              //          MyProcessorNumber, j, i, FromRegion[j].Processor,
              //          ToRegion[i].Processor);
 
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		cSndRcv[oSet+n].Sends[sends].StartIndex[dim] = LeftIndex[dim] -
		  FromRegion[j].StartIndex[dim];
		cSndRcv[oSet+n].Sends[sends].RegionDim[dim] = RightIndex[dim] -
		  LeftIndex[dim] + 1;
	      }
	      SendSize += size;
	      cSndRcv[oSet+n].Sends[sends++].Processor = j;
	    }
 
	    if (MyProcessorNumber == ToRegion[i].Processor) {
 
              //  fprintf(stderr, "CT(%"ISYM"): to: %"ISYM"->%"ISYM" from proc %"ISYM"\n",
              //          MyProcessorNumber, j, i, FromRegion[j].Processor);
 
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		cSndRcv[oSet+n].Receives[receives].StartIndex[dim] = LeftIndex[dim] -
		  ToRegion[i].StartIndex[dim];
		cSndRcv[oSet+n].Receives[receives].RegionDim[dim] = RightIndex[dim] -
		  LeftIndex[dim] + 1;
	      }
 
	      ReceiveSize += size;
	      cSndRcv[oSet+n].Receives[receives++].Processor = i;
 
	    }
 
	  } // end: if (size > 0)
 
        } // end: if (proc jump == n)
 
	cSndRcv[oSet+n].SendSize = SendSize;
	cSndRcv[oSet+n].ReceiveSize = ReceiveSize;
	cSndRcv[oSet+n].sends = sends;
	cSndRcv[oSet+n].receives = receives;
 
  } // end loop over NumberOfProcessors
  } // end if First_Pass < 6
 
  // This completes the collection and storage of the communication pattern
  // after the first 6 calls on the first pass through the FFT
 
  //
  //
 
  // Use the stored patterns in the first and all subsequent calls
 
  for (n = 0; n < NumberOfProcessors; n++) {
 
    SendSize = cSndRcv[oSet+n].SendSize;
    ReceiveSize = cSndRcv[oSet+n].ReceiveSize;
    sends = cSndRcv[oSet+n].sends;
    receives = cSndRcv[oSet+n].receives;
    Sends = cSndRcv[oSet+n].Sends;
    Receives = cSndRcv[oSet+n].Receives;
 
    /* Allocate buffer and copy data into buffer. */
 
    ReceiveBuffer = NULL;
    SendBuffer = new float[SendSize];
 
    index = 0;
 
//  fprintf(stderr, "CT(%"ISYM"): sends = %"ISYM"  SendSize = %"ISYM" receives = %"ISYM"  RecvSzs = %"ISYM"\n", MyProcessorNumber,
//	    sends, SendSize, receives, ReceiveSize);
 
    for (i = 0; i < sends; i++) {
 
      j = Sends[i].Processor;
 
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
 
      //  fprintf(stderr, "CT(%"ISYM"): MPI SS/RS = %"ISYM"/%"ISYM" From/To = %"ISYM" %"ISYM"\n",
      //  MyProcessorNumber, SendSize, ReceiveSize, FromProc, ToProc);

      Count = SendSize;
      RecvCount = ReceiveSize;
      Source = FromProc;
      Dest = ToProc;
 
      if (MPI_Sendrecv((void*) SendBuffer, Count, DataType, Dest,
	       MPI_TRANSPOSE_TAG, (void*) ReceiveBuffer, RecvCount,
	       DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
	       &status) != MPI_SUCCESS) {
	ENZO_VFAIL("Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber,
		status.MPI_ERROR)
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
 
      if (ToRegion[j].Data == NULL)
	ToRegion[j].Data = new float[ToRegion[j].RegionDim[0]*
	      ToRegion[j].RegionDim[1]*ToRegion[j].RegionDim[2]];
 
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
 
    //  fprintf(stderr, "CT(%"ISYM"): end jump %"ISYM"\n", MyProcessorNumber, n);
 
    delete [] SendBuffer;
    delete [] ReceiveBuffer;
 
  } // end: loop over processors jumps
 
 
  /* Delete FromRegion data. */
 
  //  fprintf(stderr, "CT(%"ISYM"): Deleting FromRegions\n", MyProcessorNumber);
 
  for (i = 0; i < NumberOfFromRegions; i++) {
    delete [] FromRegion[i].Data;
    FromRegion[i].Data = NULL;
  }
 
  First_Pass++;
 
  return SUCCESS;
}
 
#endif
