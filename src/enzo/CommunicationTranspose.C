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
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include "EnzoTiming.h"
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

int CommunicationBarrier(void);
#ifdef USE_MPI
int CommunicationBufferPurge(void);
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif

int TransposeRegionOverlap(region *FromRegion, region *ToRegion, int i, int j,
			   region *Sends, region *Receives, 
			   int &nrecv, int &nsend,
			   int &SendSize, int &ReceiveSize);

int NonUnigridCommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder);
int OptimizedUnigridCommunicationTranspose(region *FromRegion, 
			   int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder);
int NonBlockingCommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
		           int TransposeOrder);

#define DEBUG_NONBLOCKCT_OFF
commSndRcv *cSndRcv; /* Used in optimized transpose */
 
int CommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{
    TIMER_START("CommunicationTranspose");
    int retval;
    switch (UnigridTranspose) {
    case 0:
      retval = NonUnigridCommunicationTranspose
	(FromRegion, NumberOfFromRegions,
	 ToRegion, NumberOfToRegions, TransposeOrder);
      break;
    case 1:
      retval = OptimizedUnigridCommunicationTranspose
	(FromRegion, NumberOfFromRegions,
	 ToRegion, NumberOfToRegions, TransposeOrder);
      break;
    case 2:
      retval = NonBlockingCommunicationTranspose
	(FromRegion, NumberOfFromRegions,
	 ToRegion, NumberOfToRegions, TransposeOrder);
      break;
    default:
      ENZO_VFAIL("Invalid value for UnigridTranspose = %d", UnigridTranspose);
    } // ENDSWITCH
    TIMER_STOP("CommunicationTranspose");
    return retval;
}

int NonUnigridCommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{

  /* Declarations. */
 
  int dim, n, i, j, size, index, Zero[] = {0,0,0};
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];
  float *ReceiveBuffer, *SendBuffer;
 
  //  fprintf(stderr, "CT(%"ISYM"): start From=%"ISYM"  To=%"ISYM"\n", MyProcessorNumber,
  //	  NumberOfFromRegions, NumberOfToRegions);

  int NumberOfRegions = max(NumberOfFromRegions, NumberOfToRegions);
  region *Sends = new region[NumberOfRegions];
  region *Receives = new region[NumberOfRegions];

  //  if (NumberOfProcessors == 1) return SUCCESS;
 
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

	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, Sends, Receives,
				 receives, sends, SendSize, ReceiveSize);
 
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
 
      MPI_Request RequestHandle;
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
       
//      if (MPI_Sendrecv((void*) SendBuffer, Count, DataType, Dest,
//	       MPI_TRANSPOSE_TAG, (void*) ReceiveBuffer, RecvCount,
//	       DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
//	       &status) != MPI_SUCCESS) {
//	ENZO_VFAIL("Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber,
//		status.MPI_ERROR)
//      }

      MPI_Irecv((void*) ReceiveBuffer, RecvCount, DataType, Source, 
		MPI_TRANSPOSE_TAG, MPI_COMM_WORLD, &RequestHandle);
      MPI_Send((void*) SendBuffer, Count, DataType, Dest, 
	       MPI_TRANSPOSE_TAG, MPI_COMM_WORLD);
      MPI_Wait(&RequestHandle, &status);
 
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

      if (ToRegion[j].Data == NULL) {
	int ijk, rsize = ToRegion[j].RegionDim[0]*ToRegion[j].RegionDim[1]*ToRegion[j].RegionDim[2];
	ToRegion[j].Data = new float[rsize];
	for (ijk = 0; ijk < rsize; ijk++)
	  ToRegion[j].Data[ijk] = 0;
      }

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
};

// #endif 


int OptimizedUnigridCommunicationTranspose(
               region *FromRegion, int NumberOfFromRegions,
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
 
  // The first time this routine is called in each of 6 axes,
  // collect the static routing information for unigrid and
  // create a static table of routing information 
 
  if(First_Pass == 0)
    cSndRcv  = new commSndRcv[6*NumberOfProcessors];	



 
  if(First_Pass < 6)  {

//    fprintf(stderr, "FFT Initialization Pass %"ISYM" on proc %"ISYM" : #FromRegions =  %"ISYM"  #ToRegions =  %"ISYM"\n",
//	    First_Pass, MyProcessorNumber, NumberOfFromRegions, NumberOfToRegions);

  // Loop over processor jumps (number of processors ahead to send)

  int MaxRegions = max(NumberOfFromRegions, NumberOfToRegions);
 
  int *jtrue = new int[NumberOfFromRegions];
  int *itrue = new int[NumberOfToRegions];

  int *jfalse = new int[NumberOfFromRegions];
  int *ifalse = new int[NumberOfToRegions];

  for (n = 0; n < NumberOfProcessors; n++) {
 
    cSndRcv[oSet+n].Sends = new region[NumberOfFromRegions];
    cSndRcv[oSet+n].Receives = new region[NumberOfToRegions];

    // Initialization may be unnecessary

/*
    cSndRcv[oSet+n].sends = 0;
    cSndRcv[oSet+n].receives = 0;
    cSndRcv[oSet+n].SendSize = 0;
    cSndRcv[oSet+n].ReceiveSize = 0; 

    for (j = 0; j < NumberOfFromRegions; j++) {
      cSndRcv[oSet+n].Sends[j].Processor = 0;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
        cSndRcv[oSet+n].Sends[j].StartIndex[dim] = 0;
        cSndRcv[oSet+n].Sends[j].RegionDim[dim] = 0;
      }
    }
    for (i = 0; i < NumberOfToRegions; i++) {
      cSndRcv[oSet+n].Receives[i].Processor = 0;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
        cSndRcv[oSet+n].Receives[i].StartIndex[dim] = 0;
        cSndRcv[oSet+n].Receives[i].RegionDim[dim] = 0;
      }
    }
*/

    /* Copy regions into communication buffer (or just set buffer
       if there is only one FromRegion per processor) */
 
    sends = 0, receives = 0;
    SendSize = 0, ReceiveSize = 0;

    int jct = 0;
    int jcf = 0;

    for (j = 0; j < NumberOfFromRegions; j++) {
      if ( MyProcessorNumber == FromRegion[j].Processor ) {
        jtrue[jct] = j;
        jct = jct + 1;
      } else {
        jfalse[jcf] = j;
        jcf = jcf + 1;
      }
    }

    int ict = 0;
    int icf = 0;
    for (i = 0; i < NumberOfToRegions; i++) {
      if ( MyProcessorNumber ==   ToRegion[i].Processor ) {
        itrue[ict] = i;
        ict = ict + 1;
      } else {
        ifalse[icf] = i;
        icf = icf + 1;
      }
    }

    //    fprintf(stderr, "  Counter indices on proc %"ISYM" for n = %"ISYM": ict = %"ISYM", jct = %"ISYM"\n", MyProcessorNumber, n, ict, jct);


/*
    for (j = 0; j < NumberOfFromRegions; j++) {
      for (i = 0; i < NumberOfToRegions; i++) {

        if ((ToRegion[i].Processor - FromRegion[j].Processor +
             NumberOfProcessors) % NumberOfProcessors == n &&
            (MyProcessorNumber == FromRegion[j].Processor ||
             MyProcessorNumber ==   ToRegion[i].Processor)) {
*/

    int ii, ji;

    for ( ji = 0; ji < jct; ji++ ) {
      j = jtrue[ji];
      for ( ii = 0; ii < ict; ii++ ) {
        i = itrue[ii];

        if ( (ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) % 
	     NumberOfProcessors == n ) {

	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, 
				 cSndRcv[oSet+n].Sends, 
				 cSndRcv[oSet+n].Receives,
				 receives, sends, SendSize, ReceiveSize);
   
        } // end: if (proc jump == n)

      } // end: for ii

    } // end: for ji




    for ( ji = 0; ji < jct; ji++ ) {
      j = jtrue[ji];
      for ( ii = 0; ii < icf; ii++ ) {
        i = ifalse[ii];

        if ( (ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) % 
	     NumberOfProcessors == n ) {

	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, 
				 cSndRcv[oSet+n].Sends, 
				 cSndRcv[oSet+n].Receives,
				 receives, sends, SendSize, ReceiveSize);
   
        } // end: if (proc jump == n)

      } // end: for ii

    } // end: for ji




    for ( ji = 0; ji < jcf; ji++ ) { 
      j = jfalse[ji];
      for ( ii = 0; ii < ict; ii++ ) {
        i = itrue[ii];

        if ( (ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) % 
	     NumberOfProcessors == n ) {

	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, 
				 cSndRcv[oSet+n].Sends, 
				 cSndRcv[oSet+n].Receives,
				 receives, sends, SendSize, ReceiveSize);
   
        } // end: if (proc jump == n)

      } // end: for ii

    } // end: for ji




/*
    for ( ji = 0; ji < jcf; ji++ ) {               
      j = jfalse[ji];
      for ( ii = 0; ii < icf; ii++ ) {
        i = ifalse[ii];
      }
    }
*/




    cSndRcv[oSet+n].SendSize = SendSize;
    cSndRcv[oSet+n].ReceiveSize = ReceiveSize;
    cSndRcv[oSet+n].sends = sends;
    cSndRcv[oSet+n].receives = receives;
 
  } // end loop over NumberOfProcessors

  delete [] itrue;
  delete [] jtrue;
  delete [] ifalse;
  delete [] jfalse;

  } // end if First_Pass < 6
 
  // This completes the collection and storage of the communication pattern
  // after the first 6 calls on the first pass through the FFT




  int m, l, k;

  if ( MyProcessorNumber == -4 )
  if ( First_Pass == 5 ) {

    for ( n = 0; n < 6; n++ ) {

    for ( m = 0; m < NumberOfProcessors; m++ ) {
      fprintf(stderr, "  %"ISYM"  %"ISYM"  %"ISYM"  %"ISYM"\n",
              n, m, cSndRcv[n*NumberOfProcessors+m].ReceiveSize, cSndRcv[n*NumberOfProcessors+m].receives);
    }

    for ( m = 0; m < NumberOfProcessors; m++ ) {
      for ( l = 0; l < NumberOfToRegions; l++ ) {
        fprintf(stderr, "  %"ISYM, cSndRcv[n*NumberOfProcessors+m].Receives[l].Processor);
      }
      fprintf(stderr, "\n");
    }
/*
    for ( m = 0; m < NumberOfProcessors; m++ ) {
      for ( l = 0; l < NumberOfToRegions; l++ ) {
        for ( k = 0; k < MAX_DIMENSION; k++ ) {
           fprintf(stderr, "  %"ISYM"  %"ISYM, cSndRcv[n*NumberOfProcessors+m].Receives[l].StartIndex[k], cSndRcv[n*NumberOfProcessors+m].Receives[l].RegionDim[k]);
        }
        fprintf(stderr, "\n");
      }
    }
*/
    }



/*
    for ( m = 0; m < 6*NumberOfProcessors; m++ ) {
      for ( l = 0; l < NumberOfToRegions; l++ ) {
        for ( k = 0; k < MAX_DIMENSION; k++ ) {
           fprintf(stderr, "  %"ISYM"  %"ISYM, cSndRcv[m].Receives[l].StartIndex[k], cSndRcv[m].Receives[l].RegionDim[k]);
        }
        fprintf(stderr, "  %"ISYM, cSndRcv[m].Receives[l].Processor);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  %"ISYM"\n", cSndRcv[m].ReceiveSize);
      fprintf(stderr, "  %"ISYM"\n", cSndRcv[m].receives);
    }

    for ( m = 0; m < 6*NumberOfProcessors; m++ ) {
      for ( l = 0; l < NumberOfFromRegions; l++ ) {
        for ( k = 0; k < MAX_DIMENSION; k++ ) {
          fprintf(stderr, "  %"ISYM"  %"ISYM, cSndRcv[m].Sends[l].StartIndex[k], cSndRcv[m].Sends[l].RegionDim[k]);
        }
        fprintf(stderr, "  %"ISYM, cSndRcv[m].Sends[l].Processor);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  %"ISYM"\n", cSndRcv[m].SendSize);
      fprintf(stderr, "  %"ISYM"\n", cSndRcv[m].sends);
    }
*/

  }

  // cSndRcv[0, 6*Ncpu-1].Receives[0, NumberOfToRegions-1].StartIndex[0,2]
  // cSndRcv[0, 6*Ncpu-1].Receives[0, NumberOfToRegions-1].RegionDim[0,2]
  // cSndRcv[0, 6*Ncpu-1].Receives[0, NumberOfToRegions-1].Processor

  // cSndRcv[0, 6*Ncpu-1].Sends[0, NumberOfFromRegions-1].StartIndex[0.2]
  // cSndRcv[0, 6*Ncpu-1].Sends[0, NumberOfFromRegions-1].RegionDim[0.2]
  // cSndRcv[0, 6*Ncpu-1].Sends[0, NumberOfFromRegions-1].Processor
 
  // cSndRcv[0, 6*Ncpu-1].SendSize
  // cSndRcv[0, 6*Ncpu-1].ReceiveSize
  // cSndRcv[0, 6*Ncpu-1].sends
  // cSndRcv[0, 6*Ncpu-1].receives



 
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

      MPI_Request RequestHandle; 
      MPI_Status Status;
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

/* 
      if (MPI_Sendrecv((void*) SendBuffer, Count, DataType, Dest,
	       MPI_TRANSPOSE_TAG, (void*) ReceiveBuffer, RecvCount,
	       DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
	       &Status) != MPI_SUCCESS) {
	ENZO_VFAIL("Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber,
		Status.MPI_ERROR)
      }
*/

      MPI_Irecv((void*) ReceiveBuffer, RecvCount, DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD, &RequestHandle);
      MPI_Send((void*) SendBuffer, Count, DataType, Dest, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD);
      MPI_Wait(&RequestHandle, &Status);

 
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
 
    //fprintf(stderr, "CT(%"ISYM"): receives = %"ISYM"\n", MyProcessorNumber, receives);
 
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
 
    //fprintf(stderr, "CT(%"ISYM"): end jump %"ISYM"\n", MyProcessorNumber, n);
 
    delete [] SendBuffer;
    delete [] ReceiveBuffer;
 
  } // end: loop over processors jumps
 
 
  /* Delete FromRegion data. */
 
  //fprintf(stderr, "CT(%"ISYM"): Deleting FromRegions\n", MyProcessorNumber);
 
  for (i = 0; i < NumberOfFromRegions; i++) {
    delete [] FromRegion[i].Data;
    FromRegion[i].Data = NULL;
  }
 
  First_Pass++;
 
  return SUCCESS;

}

/* Non-blocking version of CommunicationTranspose where we process the
   received data every PROCS_PER_LOOP cycles. */ 

#define PROCS_PER_LOOP 128
int NonBlockingCommunicationTranspose(region *FromRegion, int NumberOfFromRegions,
			   region *ToRegion, int NumberOfToRegions,
			   int TransposeOrder)
{
#ifdef USE_MPI
  /* Declarations. */
 
  int dim, n, ni, i, ii, j, jj, size, index, Zero[] = {0,0,0};
  int LeftIndex[MAX_DIMENSION], RightIndex[MAX_DIMENSION];
  float *ReceiveBuffer[PROCS_PER_LOOP], *SendBuffer[PROCS_PER_LOOP];
  bool ReceiveMode;
  int sends, receives, request;
  int SendSize, ReceiveSize;
  int NumberOfRequests = 0;

  MPI_Request RequestHandle[PROCS_PER_LOOP];
  MPI_Status ListOfStatuses[PROCS_PER_LOOP];
  MPI_Arg ListOfIndices[PROCS_PER_LOOP];
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg RecvCount;
  MPI_Arg Source;
  MPI_Arg Dest;
  MPI_Arg error_code;
  char error_string[1024];
  MPI_Arg length_of_error_string, error_class;
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
 
#ifdef DEBUG_NONBLOCKCT
    fprintf(stderr, "CT(%"ISYM"): start From=%"ISYM"  To=%"ISYM"\n", 
	    MyProcessorNumber, NumberOfFromRegions, NumberOfToRegions);
#endif

  int NumberOfRegions = max(NumberOfFromRegions, NumberOfToRegions);
  commSndRcv commNB[PROCS_PER_LOOP];
  region *Sends, *Receives;

  for (i = 0; i < PROCS_PER_LOOP; i++) {
    commNB[i].Sends = new region[NumberOfRegions];
    commNB[i].Receives = new region[NumberOfRegions];
    //MPI_Request_free(RequestHandle+i);
  }

  //  if (NumberOfProcessors == 1) return SUCCESS;

  /* Store which From/To Regions this processor will process */

  int jct = 0, jcf = 0, ict = 0, icf = 0;
  int *jtrue = new int[NumberOfFromRegions];
  int *itrue = new int[NumberOfToRegions];
  int *jfalse = new int[NumberOfFromRegions];
  int *ifalse = new int[NumberOfToRegions];

  for (j = 0; j < NumberOfFromRegions; j++)
    if (MyProcessorNumber == FromRegion[j].Processor)
      jtrue[jct++] = j;
    else
      jfalse[jcf++] = j;
  for (i = 0; i < NumberOfToRegions; i++)
    if (MyProcessorNumber == ToRegion[i].Processor)
      itrue[ict++] = i;
    else
      ifalse[icf++] = i;
 
  /* Loop over processor jumps (number of processors ahead to send). */
 
  for (n = 0; n < NumberOfProcessors; n++) {
 
    /* Receive every PROCS_PER_LOOP (and last) cycle.  The first cycle
       is special because there is no communication. */

    ni = (max(n-1,0)) % PROCS_PER_LOOP;
    ReceiveMode = (ni == PROCS_PER_LOOP-1 || n == NumberOfProcessors-1 || n==0);

    sends = receives = 0;
    SendSize = ReceiveSize = 0;
    Sends = commNB[ni].Sends;
    Receives = commNB[ni].Receives;

    /* Copy regions into communication buffer (or just set buffer if
       there is only one FromRegion per processor).  Check only
       regions that involve this processor. */

//    for (j = 0; j < NumberOfFromRegions; j++)
//      for (i = 0; i < NumberOfToRegions; i++)
//        if ((ToRegion[i].Processor - FromRegion[j].Processor +
//             NumberOfProcessors) % NumberOfProcessors == n &&
//	    (MyProcessorNumber == FromRegion[j].Processor ||
//	     MyProcessorNumber ==   ToRegion[i].Processor)) {

    for (jj = 0; jj < jct; jj++) {
      j = jtrue[jj];
      for (ii = 0; ii < ict; ii++) {
	i = itrue[ii];
	if ((ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) %
	    NumberOfProcessors == n)
	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, Sends, Receives,
				 receives, sends, SendSize, ReceiveSize);
      } // ENDFOR ii
    } // ENDFOR jj

    for (jj = 0; jj < jct; jj++) {
      j = jtrue[jj];
      for (ii = 0; ii < icf; ii++) {
	i = ifalse[ii];
	if ((ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) %
	    NumberOfProcessors == n)
	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, Sends, Receives,
				 receives, sends, SendSize, ReceiveSize);
      } // ENDFOR ii
    } // ENDFOR jj

    for (jj = 0; jj < jcf; jj++) {
      j = jfalse[jj];
      for (ii = 0; ii < ict; ii++) {
	i = itrue[ii];
	if ((ToRegion[i].Processor - FromRegion[j].Processor + NumberOfProcessors) %
	    NumberOfProcessors == n)
	  TransposeRegionOverlap(FromRegion, ToRegion, i, j, Sends, Receives,
				 receives, sends, SendSize, ReceiveSize);
      } // ENDFOR ii
    } // ENDFOR jj

    /* Store data into structure for the last PROCS_PER_LOOP cycles */

    commNB[ni].sends = sends;
    commNB[ni].receives = receives;
    commNB[ni].SendSize = SendSize;
    commNB[ni].ReceiveSize = ReceiveSize;
 
    /* Allocate buffer and copy data into buffer. */
 
    ReceiveBuffer[ni] = NULL;
    SendBuffer[ni] = new float[SendSize];
 
    index = 0;
 
#ifdef DEBUG_NONBLOCKCT
    fprintf(stderr, "CT(%"ISYM"): sends = %"ISYM"  SendSize = %"ISYM"\n", 
	    MyProcessorNumber, sends, SendSize);
#endif
 
    for (i = 0; i < sends; i++) {
      j = Sends[i].Processor;
      if (FromRegion[j].Data == NULL) {
	int ijk, rsize;
	rsize = FromRegion[j].RegionDim[0] * 
	  FromRegion[j].RegionDim[1] * 
	  FromRegion[j].RegionDim[2];

	// allocate FROM field
	FromRegion[j].Data = new float[rsize];
	for (ijk = 0; ijk < rsize; ijk++)
	  FromRegion[j].Data[ijk] = 0;
      }
      if (TransposeOrder == TRANSPOSE_REVERSE)
	FORTRAN_NAME(copy3drt)(FromRegion[j].Data, SendBuffer[ni]+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1,
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1,
			   Sends[i].StartIndex+2);
      else
	FORTRAN_NAME(copy3d)(FromRegion[j].Data, SendBuffer[ni]+index,
		       FromRegion[j].RegionDim, FromRegion[j].RegionDim+1,
			   FromRegion[j].RegionDim+2,
		       Sends[i].RegionDim, Sends[i].RegionDim+1,
			   Sends[i].RegionDim+2,
		       Zero, Zero+1, Zero+2,
		       Sends[i].StartIndex, Sends[i].StartIndex+1,
			   Sends[i].StartIndex+2);
      index += Sends[i].RegionDim[0] * Sends[i].RegionDim[1] * 
	Sends[i].RegionDim[2];
    }
 
    /* shift buffer by n processors. */
 
    if (n > 0) {
 
      ReceiveBuffer[ni] = new float[ReceiveSize];
 
      int ToProc = (MyProcessorNumber + n) % NumberOfProcessors;
      int FromProc = (MyProcessorNumber - n + NumberOfProcessors) %
	NumberOfProcessors;
 
#ifdef MPI_INSTRUMENTATION
      starttime = MPI_Wtime();
#endif
 
#ifdef DEBUG_NONBLOCKCT
      fprintf(stderr, "CT(%"ISYM"): MPI SS/RS = %"ISYM"/%"ISYM" From/To = %"ISYM" %"ISYM"\n",
	      MyProcessorNumber, SendSize, ReceiveSize, FromProc, ToProc);
#endif

      Count = SendSize;
      RecvCount = ReceiveSize;
      Source = FromProc;
      Dest = ToProc;
       
//      if (MPI_Sendrecv((void*) SendBuffer, Count, DataType, Dest,
//	       MPI_TRANSPOSE_TAG, (void*) ReceiveBuffer, RecvCount,
//	       DataType, Source, MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
//	       &status) != MPI_SUCCESS) {
//	ENZO_VFAIL("Proc %"ISYM" MPI_Sendrecv error %"ISYM"\n", MyProcessorNumber,
//		status.MPI_ERROR)
//      }

      /* Post receive call */

      MPI_Irecv((void*) ReceiveBuffer[ni], RecvCount, DataType, Source, 
		MPI_TRANSPOSE_TAG, MPI_COMM_WORLD, RequestHandle+ni);

      /* Buffered send of the data */

      CommunicationBufferedSend(SendBuffer[ni], Count, DataType, Dest,
				MPI_TRANSPOSE_TAG, MPI_COMM_WORLD,
				BUFFER_IN_PLACE);

//      MPI_Send((void*) SendBuffer[ni], Count, DataType, Dest, 
//	       MPI_TRANSPOSE_TAG, MPI_COMM_WORLD);
//      MPI_Wait(&RequestHandle, &status);
 
#ifdef MPI_INSTRUMENTATION
      endtime = MPI_Wtime();
      timer[14] += endtime-starttime;
      counter[14] ++;
      CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
    } else {
      //MPI_Request_free(RequestHandle+ni);
      ReceiveBuffer[ni] = SendBuffer[ni];
      SendBuffer[ni] = NULL;
    }

    // Every cycle has a request (valid or free)
    NumberOfRequests++;
 
    /* Copy from communication buffer back to regions if we're
       processing receives this cycle. */
 
    MPI_Arg TotalCompletedRequests, CompletedRequests;
    if (ReceiveMode) {

      TotalCompletedRequests = 0;

#ifdef DEBUG_NONBLOCKCT
      printf("CT(%"ISYM"): TotalCompletedRequests = %d, NumberOfRequests = %d\n",
	     MyProcessorNumber, TotalCompletedRequests, 
	     NumberOfRequests);
#endif
      
      while (TotalCompletedRequests < NumberOfRequests) {

	if (n > 0) {
	  CompletedRequests = 0;
	  MPI_Waitsome(NumberOfRequests, RequestHandle, &CompletedRequests,
		       ListOfIndices, ListOfStatuses);
#ifdef DEBUG_NONBLOCKCT
	  printf("P%d: CompletedRequests = %d/%d\n", MyProcessorNumber, 
		 CompletedRequests, NumberOfRequests);
#endif

	  /* Check for any errors */

#ifdef UNUSED
	  for (i = 0; i < NumberOfRequests; i++) {
	    index = ListOfIndices[i];
	    if (ListOfStatuses[index].MPI_ERROR != 0) {
	      fprintf(stderr, "CommunicationTranspose: MPI Error on request %d, "
		      "n = %d, error = %d\n", i, n, ListOfStatuses[i].MPI_ERROR);
	      MPI_Error_class(ListOfStatuses[index].MPI_ERROR, &error_class);
	      MPI_Error_string(error_class, error_string, &length_of_error_string);
	      fprintf(stderr, "P%d: %s\n", MyProcessorNumber, error_string);
	      MPI_Error_string(ListOfStatuses[index].MPI_ERROR, error_string, 
			       &length_of_error_string);
	      fprintf(stderr, "P%d: %s\n", MyProcessorNumber, error_string);
	      ENZO_FAIL("MPI Error in CommunicationTranspose!\n");
	    }
	  }
#endif

	} // ENDIF n > 0

	for (request = 0; request < NumberOfRequests; request++) {

	  /* First cycle only has one call.  Break if we've already
	     processed the first call, otherwise always process it. */

	  if (n == 0 && request > 0) break;
	  if ((RequestHandle[request] == MPI_REQUEST_NULL &&
	      ReceiveBuffer[request] != NULL) || 
	      (n == 0 && request == 0)) {

#ifdef DEBUG_NONBLOCKCT
	    fprintf(stderr, "CT(%"ISYM"): request = %d, receives = %"ISYM"\n", 
		    MyProcessorNumber, request, receives);
#endif

	    index = 0;
	    Receives = commNB[request].Receives;
	    receives = commNB[request].receives;

	    for (i = 0; i < receives; i++) {

	      j = Receives[i].Processor;

	      if (ToRegion[j].Data == NULL) {
		int ijk, rsize;
		rsize = ToRegion[j].RegionDim[0] * 
		  ToRegion[j].RegionDim[1] * 
		  ToRegion[j].RegionDim[2];
		ToRegion[j].Data = new float[rsize];
		for (ijk = 0; ijk < rsize; ijk++)
		  ToRegion[j].Data[ijk] = 0;
	      }

	      if (TransposeOrder == TRANSPOSE_FORWARD)

		FORTRAN_NAME(copy3dft)(ReceiveBuffer[request]+index, 
				 ToRegion[j].Data,
				 Receives[i].RegionDim, Receives[i].RegionDim+1,
				 Receives[i].RegionDim+2,
				 ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
				 ToRegion[j].RegionDim+2,
				 Receives[i].StartIndex, Receives[i].StartIndex+1,
				 Receives[i].StartIndex+2,
				 Zero, Zero+1, Zero+2);
	      else
		FORTRAN_NAME(copy3d)(ReceiveBuffer[request]+index, 
			       ToRegion[j].Data,
			       Receives[i].RegionDim, Receives[i].RegionDim+1,
			       Receives[i].RegionDim+2,
			       ToRegion[j].RegionDim, ToRegion[j].RegionDim+1,
			       ToRegion[j].RegionDim+2,
			       Receives[i].StartIndex, Receives[i].StartIndex+1,
			       Receives[i].StartIndex+2,
			       Zero, Zero+1, Zero+2);
	      index += Receives[i].RegionDim[0]*Receives[i].RegionDim[1]*
		Receives[i].RegionDim[2];
	    } // ENDFOR receives

	    delete [] ReceiveBuffer[request];
	    ReceiveBuffer[request] = NULL;
	    TotalCompletedRequests++;

	  } // ENDIF completed request
	} // ENDFOR requests
#ifdef DEBUG_NONBLOCKCT
	printf("CT(%"ISYM"): (n=%d) -- completed %d out of %d requests\n",
	       MyProcessorNumber, n, TotalCompletedRequests, NumberOfRequests);
#endif
      } // ENDWHILE

      /* Clean up (SendBuffer will be deleted in CommunicationBufferedSend). */
 
//      delete [] SendBuffer;
//      for (i = 0; i < PROCS_PER_LOOP; i++)
//	MPI_Request_free(RequestHandle+i);
//	delete [] ReceiveBuffer[i];

      NumberOfRequests = 0;

    } // ENDIF ReceiveMode

#ifdef DEBUG_NONBLOCKCT
    fprintf(stderr, "CT(%"ISYM"): end jump %"ISYM"\n", MyProcessorNumber, n);
#endif
 
  } // end: loop over processors jumps
 
  /* Delete FromRegion data. */
 
  for (i = 0; i < NumberOfFromRegions; i++) {
    delete [] FromRegion[i].Data;
    FromRegion[i].Data = NULL;
  }
 
  /* Clean up. */

  delete [] itrue;
  delete [] jtrue;
  delete [] ifalse;
  delete [] jfalse;

  for (i = 0; i < PROCS_PER_LOOP; i++) {
    delete [] commNB[i].Sends;
    delete [] commNB[i].Receives;
  }

  CommunicationBarrier();
  CommunicationBufferPurge();

#else
  ENZO_FAIL("UnigridTranspose = 2 can only be used with use-mpi-yes.");
#endif /* USE_MPI */
 
  return SUCCESS;
};
