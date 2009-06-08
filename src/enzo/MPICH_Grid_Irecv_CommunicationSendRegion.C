/***********************************************************************
/
/  GRID CLASS (SEND FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
// function prototypes
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
 
int grid::CommunicationSendRegion(grid *ToGrid, int ToProcessor,int SendField,
			      int NewOrOld, int RegionStart[], int RegionDim[])
{
 
  MPI_Request  RequestHandle;
  MPI_Status Status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg Source;

 
  int index, field, dim, Zero[] = {0, 0, 0};
 
  // Compute size of region to transfer
 
  int NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
                       ((NewOrOld == NEW_AND_OLD)? 2 : 1);
 
  if (SendField == ACCELERATION_FIELDS)
    NumberOfFields = GridRank;
 
  int RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
  int TransferSize = RegionSize * NumberOfFields;
 
  // Allocate buffer
 
  float *buffer = NULL;
 
  if (MyProcessorNumber == ProcessorNumber || MyProcessorNumber == ToProcessor)
    buffer = new float[TransferSize];
 
  // If this is the from processor, pack fields
 
  if (MyProcessorNumber == ProcessorNumber) {
 
//  printf("SendRegion: RegionStart = %"ISYM" %"ISYM" %"ISYM"\n", RegionStart[0], RegionStart[1], RegionStart[2]);
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(BaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(OldBaryonField[field], &buffer[index],
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       RegionStart, RegionStart+1, RegionStart+2);
	  index += RegionSize;
      }
 
    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES)
      FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, buffer,
			   GravitatingMassFieldParticlesDimension,
			   GravitatingMassFieldParticlesDimension+1,
			   GravitatingMassFieldParticlesDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == GRAVITATING_MASS_FIELD)
      FORTRAN_NAME(copy3d)(GravitatingMassField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == POTENTIAL_FIELD)
      FORTRAN_NAME(copy3d)(PotentialField, buffer,
			   GravitatingMassFieldDimension,
			   GravitatingMassFieldDimension+1,
			   GravitatingMassFieldDimension+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   RegionStart, RegionStart+1, RegionStart+2);
 
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	FORTRAN_NAME(copy3d)(AccelerationField[dim], &buffer[index],
			     GridDimension, GridDimension+1, GridDimension+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     RegionStart, RegionStart+1, RegionStart+2);
	index += RegionSize;
      }
  }
 
  // Send buffer
 
#ifdef USE_MPI
 
  // Only send if processor numbers are not identical
 
  if (ProcessorNumber != ToProcessor) {
 
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
 
    if (MyProcessorNumber == ProcessorNumber) {
#ifdef MPI_INSTRUMENTATION
      if (traceMPI) fprintf(tracePtr, "CSR Sending %"ISYM" floats from %"ISYM" to %"ISYM"\n", TransferSize,
	      MyProcessorNumber, ToProcessor);
#endif
      CommunicationBufferedSend(buffer, TransferSize, DataType, ToProcessor,
	       MPI_SENDREGION_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    }
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[5] += endtime-starttime;
    counter[5] ++;
    timer[6] += double(TransferSize);
    RecvComm += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
  }
 
  // If this is the to processor then post the receive
  // Only post if processor numbers are not identical
 
  if (ProcessorNumber != ToProcessor) {
 
    if (MyProcessorNumber == ToProcessor) {
#ifdef MPI_INSTRUMENTATION
      if (traceMPI) fprintf(tracePtr, "CSR Post Irecv %"ISYM" %"ISYM" %"ISYM" %"ISYM" = %"ISYM"\n", RegionDim[0],RegionDim[1],RegionDim[2],NumberOfFields,TransferSize);
      if (traceMPI) fprintf(tracePtr, "CSR Post Irecv for %"ISYM" floats at %"ISYM" from %"ISYM"\n", TransferSize,
              MyProcessorNumber, ProcessorNumber);
#endif
 
#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif

      Count = TransferSize;
      Source = ProcessorNumber;
 
      MPI_Recv(buffer, Count, DataType, Source, MPI_SENDREGION_TAG, MPI_COMM_WORLD, &Status);
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[5] += endtime-starttime;
    RecvComm += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
    }
 
  }
 
#endif /* USE_MPI */
 
  /* If this is the to processor, unpack fields. */
 
  if (MyProcessorNumber == ToProcessor) {
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->BaryonField[field];
	  ToGrid->BaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < max(NumberOfBaryonFields, SendField+1); field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  delete ToGrid->OldBaryonField[field];
	  ToGrid->OldBaryonField[field] = new float[RegionSize];
	  FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       Zero, Zero+1, Zero+2);
	  index += RegionSize;
	}
 
    if (SendField == GRAVITATING_MASS_FIELD_PARTICLES) {
      delete ToGrid->GravitatingMassFieldParticles;
      ToGrid->GravitatingMassFieldParticles = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassFieldParticles,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == GRAVITATING_MASS_FIELD) {
      delete ToGrid->GravitatingMassField;
      ToGrid->GravitatingMassField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->GravitatingMassField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == POTENTIAL_FIELD) {
      delete ToGrid->PotentialField;
      ToGrid->PotentialField = new float[RegionSize];
      FORTRAN_NAME(copy3d)(buffer, ToGrid->PotentialField,
			   RegionDim, RegionDim+1, RegionDim+2,
			   RegionDim, RegionDim+1, RegionDim+2,
			   Zero, Zero+1, Zero+2,
			   Zero, Zero+1, Zero+2);
    }
 
    if (SendField == ACCELERATION_FIELDS)
      for (dim = 0; dim < GridRank; dim++) {
	delete ToGrid->AccelerationField[dim];
	ToGrid->AccelerationField[dim] = new float[RegionSize];
	FORTRAN_NAME(copy3d)(&buffer[index], ToGrid->AccelerationField[dim],
			     RegionDim, RegionDim+1, RegionDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     Zero, Zero+1, Zero+2);
	index += RegionSize;
      }
			
  }
 
  if (MyProcessorNumber == ToProcessor)
    delete [] buffer;
 
  return SUCCESS;
}
