/***********************************************************************
/
/  GRID CLASS (SEND PARTICLES FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       March, 2004
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
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
#include "CommunicationUtilities.h"

/* function prototypes */
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
static int FirstTimeCalled = TRUE;
static MPI_Datatype ParticleDataType;
#endif /* USE_MPI */

void my_exit(int status);
 
/* Send particle from this grid to ToGrid on processor ToProcessor, using
   FromNumber particles counting from FromStart.  Place into ToGrid at
   particle number ToStart. If ToStart = -1, then add to end. */
 
int grid::CommunicationSendParticles(grid *ToGrid, int ToProcessor,
				    int FromStart, int FromNumber, int ToStart)
{
#ifdef USE_MPI 
  int i, j, dim, index;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

//  if (MyProcessorNumber != ProcessorNumber && 
//      MyProcessorNumber != ToProcessor)
//    return SUCCESS;
 
  if (FromNumber == 0)
    return SUCCESS;
 
  /* Compute size of region to transfer. */
 
//  int NumberOfFields = 1+GridRank+NumberOfParticleAttributes;
//  int RegionSize = FromNumber;
  int TransferSize = FromNumber;
 
  //  fprintf(stderr, "P(%"ISYM"): sending %"ISYM" from %"ISYM" -> %"ISYM" (%"ISYM" %"ISYM")\n",
  //  	  MyProcessorNumber, FromNumber, ProcessorNumber, ToProcessor,
  //  	  FromStart, ToStart);
 
  /* Allocate buffer. */
 
  particle_data *buffer = NULL;
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (particle_data*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
    buffer = new particle_data[TransferSize];
 
  /* If this is the from processor, pack fields. */
 
  if (MyProcessorNumber == ProcessorNumber) {

    int FromEnd = FromStart + FromNumber;

    for (dim = 0; dim < GridRank; dim++) {
      index = 0;
      for (i = FromStart; i < FromEnd; i++, index++) {
	buffer[index].pos[dim] = ParticlePosition[dim][i];
	buffer[index].vel[dim] = ParticleVelocity[dim][i];
      }
    }

    index = 0;
    for (i = FromStart; i < FromEnd; i++, index++) {
      buffer[index].mass = ParticleMass[i];
      buffer[index].id = ParticleNumber[i];
      buffer[index].type = ParticleType[i];
    } // ENDFOR particles

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      index = 0;
      for (i = FromStart; i < FromEnd; i++, index++)
	buffer[index].attribute[j] = ParticleAttribute[j][i];
    }
 
  } // end: if (MyProcessorNumber)
 
  /* Allocate Number field on from processor. */
 
  FLOAT *TempPos[MAX_DIMENSION];
  float  *TempVel[MAX_DIMENSION], *TempMass,
        *TempAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT *TempNumber;
  int NewNumber = FromNumber, *TempType;
  if (ToStart == -1)
    NewNumber += ToGrid->NumberOfParticles;
 
  if (MyProcessorNumber == ToProcessor) {

    /* If ToStart == -1, then add new particles to end (first, copy to buffer,
       then allocate and copy back. */
 
    if (ToStart == -1) {
      TempMass = ToGrid->ParticleMass;
      TempNumber = ToGrid->ParticleNumber;
      TempType = ToGrid->ParticleType;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	TempPos[dim] = ToGrid->ParticlePosition[dim];
	TempVel[dim] = ToGrid->ParticleVelocity[dim];
      }
      for (j = 0; j < NumberOfParticleAttributes; j++)
	TempAttribute[j] = ToGrid->ParticleAttribute[j];
      ToGrid->ParticleNumber = NULL;  // signal that we should reallocate
    }
 
    /* If unallocated, then allocate. */
 
    if (ToGrid->ParticleNumber == NULL) {
      ToGrid->AllocateNewParticles(NewNumber);
      if (ToStart > 0) {
	ENZO_VFAIL("Unallocated Number, yet FromStart = %"ISYM"\n", FromStart)
      }
    }
 
    /* If adding to end, then copy and delete old fields. */
 
    if (ToStart == -1) {
      for (i = 0; i < ToGrid->NumberOfParticles; i++) {
	ToGrid->ParticleNumber[i] = TempNumber[i];
	ToGrid->ParticleMass[i]   = TempMass[i];
	ToGrid->ParticleType[i]   = TempType[i];
      }
      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < ToGrid->NumberOfParticles; i++) {
	  ToGrid->ParticlePosition[dim][i] = TempPos[dim][i];
	  ToGrid->ParticleVelocity[dim][i] = TempVel[dim][i];
	}
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < ToGrid->NumberOfParticles; i++)
	  ToGrid->ParticleAttribute[j][i] = TempAttribute[j][i];
	
      delete [] TempNumber;
      delete [] TempMass;
      delete [] TempType;
      for (dim = 0; dim < GridRank; dim++) {
	delete [] TempPos[dim];
	delete [] TempVel[dim];
      }
      for (j = 0; j < NumberOfParticleAttributes; j++)
	delete [] TempAttribute[j];
      ToStart = ToGrid->NumberOfParticles;
    }
 
  } // end: if (MyProcessorNumber == ToProcessor)
 
  ToGrid->NumberOfParticles = max(NewNumber, ToGrid->NumberOfParticles);
 
  /* Send buffer. */
 
  /* only send if processor numbers are not identical */
 
  if (ProcessorNumber != ToProcessor) {
 
    MPI_Status status;
    MPI_Arg PCount, Count = TransferSize;
    MPI_Arg Source = ProcessorNumber;
    MPI_Arg Dest = ToProcessor;
    MPI_Arg stat;

    if (FirstTimeCalled) {
      PCount = sizeof(particle_data);
      //  fprintf(stderr, "Size of ParticleMoveList %"ISYM"\n", Count);
      stat = MPI_Type_contiguous(PCount, MPI_BYTE, &ParticleDataType);
      stat |= MPI_Type_commit(&ParticleDataType);
      if (stat != MPI_SUCCESS) ENZO_FAIL("");
      FirstTimeCalled = FALSE;
    }

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif

    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(buffer, Count, ParticleDataType,
				Dest, MPI_SENDPART_TAG, MPI_COMM_WORLD,
				BUFFER_IN_PLACE);

    if (MyProcessorNumber == ToProcessor) {

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, Count, ParticleDataType, Source,
		  MPI_SENDPART_TAG, MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

//	printf("Posting receive from P%"ISYM" for %"ISYM" particles in "
//	       "comm index %"ISYM"\n", Source, TransferSize, 
//	       CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 14;
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = FromStart;
	CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = FromNumber;
	CommunicationReceiveArgumentInt[2][CommunicationReceiveIndex] = ToStart;

	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      }

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	MPI_Recv(buffer, Count, ParticleDataType, Source,
		 MPI_SENDPART_TAG, MPI_COMM_WORLD, &status);

    } // ENDIF (MyProcessorNumber == ToProcessor)
 
#ifdef MPI_INSTRUMENTATION
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    RecvComm += endtime - starttime;
    CommunicationTime += endtime - starttime;
#endif /* MPI_INSTRUMENTATION */
 
  } // end: if (ProcessorNumber != ToProcessor)
 
  /* If this is the to processor, unpack fields. */
 
  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    int ToEnd = ToStart + FromNumber;
    
    index = 0;
    for (i = ToStart; i < ToEnd; i++, index++) {
      ToGrid->ParticleMass[i] = buffer[index].mass;
      ToGrid->ParticleType[i] = buffer[index].type;
      ToGrid->ParticleNumber[i] = buffer[index].id;
    }

    for (dim = 0; dim < GridRank; dim++) {
      index = 0;
      for (i = ToStart; i < ToEnd; i++, index++) {
	ToGrid->ParticlePosition[dim][i] = buffer[index].pos[dim];
	ToGrid->ParticleVelocity[dim][i] = buffer[index].vel[dim];
      }
    }

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      index = 0;
      for (i = ToStart; i < ToEnd; i++, index++)
	ToGrid->ParticleAttribute[j][i] = buffer[index].attribute[j];
    }

    /* Only delete the buffer if we're in receive mode (in send mode
       it will be deleted by CommunicationBufferedSend and if we're in
       post-receive mode then it will be deleted when we get to
       receive-mode). */

    delete [] buffer;
    			
  } // end: if (MyProcessorNumber...)

 
#endif /* USE_MPI */ 
  return SUCCESS;
}
 
