/***********************************************************************
/
/  GRID CLASS (SEND ACTIVE PARTICLES FROM REAL GRID TO 'FAKE' 
/              (REPLICATED) GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  John Wise -- re-purposed for active particles
/  date:       December, 2011
/  modified2:  Stephen Skory -- fixed several big errors
/  date:       Sept, 2012
/
/  NOTES:  Adapted from grid::CommunicationSendParticles().
/
************************************************************************/

#ifdef USE_MPI
#endif /* USE_MPI */

#include "preincludes.h"
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
#include "ActiveParticle.h"

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */


/* Send active particle from this grid to ToGrid on processor
   ToProcessor, using FromNumber particles counting from FromStart.
*/

int grid::CommunicationSendActiveParticles(
        grid *ToGrid, int ToProcessor, bool DeleteParticles)
{
#ifdef USE_MPI

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (NumberOfActiveParticles == 0)
    return SUCCESS;

  char *buffer;
  Eint32 position = 0;
  int npart, i, j, NumberToSend, type, dim, index;
  int element_size, header_size, ap_id;
  int *type_element_size, *type_count;
  int type_count_index;
  int SendNumberOfActiveParticles;
  ActiveParticleType_info *ap_info;
  ActiveParticleList<ActiveParticleType> NewParticles;
  ActiveParticleList<ActiveParticleType> SendParticles;

  MPI_Status status;
  MPI_Arg Count, Source, Dest;

  /* Serial case */

  if (NumberOfProcessors == 1) {
    ToGrid->AddActiveParticles(this->ActiveParticles, 0, this->ActiveParticles.size());
    if (ToGrid != this && DeleteParticles == true)
    {
      this->DeleteActiveParticles();
    }
    return SUCCESS;
  } // ENDIF serial case

  type_element_size = new int[EnabledActiveParticlesCount];
  SendNumberOfActiveParticles = NumberOfActiveParticles;
#ifdef USE_MPI
    if (CommunicationDirection == COMMUNICATION_RECEIVE)
	  type_count = (int*) CommunicationReceiveBuffer[CommunicationReceiveIndex++];
    else
#endif
      type_count = new int[EnabledActiveParticlesCount];

  for (type = 0; type < EnabledActiveParticlesCount; type++) {

      // Here we find out how many of which type of particle we have,
      // and also record how big they are. Each processor knows the
      // type element_size, but only the from processor knows how much
      // of each type is in this grid.

      ap_info = EnabledActiveParticles[type];

      header_size = ap_info->ReturnHeaderSize();
      element_size = ap_info->ReturnElementSize();
      type_element_size[type] = element_size;

      if (MyProcessorNumber == ProcessorNumber) {
	    NumberToSend = 0;
	    for (i = 0; i < NumberOfActiveParticles; i++) {
	      if (ActiveParticles[i]->ReturnType() == type) {
	        NumberToSend++;
	      }
        type_count[type] = NumberToSend;
        }
      } // myproc = proc

  } // end for particle type

  // Just once we set up a non-blocking send/recv for the
  // numbers of each type of AP
  if (ProcessorNumber != ToProcessor) {
    Count = EnabledActiveParticlesCount;
    Source = ProcessorNumber;
    Dest = ToProcessor;
    
    // send 
    if ((MyProcessorNumber == ProcessorNumber) && 
       ((CommunicationDirection == COMMUNICATION_SEND_RECEIVE) || 
        (CommunicationDirection == COMMUNICATION_SEND))) {
      CommunicationBufferedSend(type_count, Count, IntDataType, 
        Dest, MPI_SENDAP_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    }
    
    // recv
    if (MyProcessorNumber == ToProcessor) { 
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
        MPI_Irecv(type_count, Count, IntDataType, Source, MPI_SENDAP_TAG,
          MPI_COMM_WORLD, CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
      
	    CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	    CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	    CommunicationReceiveCallType[CommunicationReceiveIndex] = 22;

	    CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) type_count;
	    CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	      CommunicationReceiveCurrentDependsOn;
	    CommunicationReceiveIndex++;
	  }
	  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
	    MPI_Recv(type_count, Count, IntDataType, Source, MPI_SENDAP_TAG,
          MPI_COMM_WORLD, &status);
      }
    } // myproc == toproc
  } // procnum != toproc

  // now we loop over all the types again.
  for (type = 0; type < EnabledActiveParticlesCount; type++) {

    ap_info = EnabledActiveParticles[type];

#ifdef USE_MPI
    if (CommunicationDirection == COMMUNICATION_RECEIVE)
	  buffer = (char*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
    else
#endif
	  buffer = new char[NumberOfActiveParticles * type_element_size[type]];
      // Above we are making the largest buffer we'll ever need, which
      // eliminates the need for the receiving processor to know how many of
      // each type of particle to make more custom-tailored buffer(s).

    /* If this is the from processor, pack fields and then delete the 
       sent active particles. */
  if ((CommunicationDirection == COMMUNICATION_SEND) || 
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)) {
    if (MyProcessorNumber == ProcessorNumber) {
      ActiveParticleList<ActiveParticleType> SendParticles(type_count[type]);
      position = 0;
      ap_id = ap_info->GetEnabledParticleID();
        for (i = 0; i < NumberOfActiveParticles; i++) {
          if (ActiveParticles[i]->ReturnType() == type) {
	        SendParticles.copy_and_insert(*ActiveParticles[i]);
	        SendParticles[position]->SetGridID(ToGrid->GetGridID());
	        position++;
	  } // for type
        } // for NumberOfAPs
      // Now we can fill the buffer because we have a list of only this type.
      ap_info->FillBuffer(SendParticles, position, buffer);
      // We can clean up immediately.
      if (DeleteParticles == true)
        this->DeleteActiveParticles();
    } // myproc == proc
  } // if sending
    
  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {
    MPI_Status status;
    Count = SendNumberOfActiveParticles * type_element_size[type];
    

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
    if (MyProcessorNumber == ProcessorNumber && 
      ((CommunicationDirection == COMMUNICATION_SEND_RECEIVE) ||
       (CommunicationDirection == COMMUNICATION_SEND))) {
      // send the actual particle data.
      CommunicationBufferedSend(buffer, Count, MPI_PACKED, 
				Dest, MPI_SENDAP_TAG + 1 + type, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);
	}

    if (MyProcessorNumber == ToProcessor) {
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      	MPI_Irecv(buffer, Count, MPI_PACKED, Source,
		  MPI_SENDAP_TAG + 1 + type, MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 22;

	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      }

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) 
	MPI_Recv(buffer, Count, MPI_PACKED, Source,
		 MPI_SENDAP_TAG + 1 + type, MPI_COMM_WORLD, &status);

    } // ENDIF MyProcessorNumber == ToProcessor

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    timer[8] += double(Count);
    timer[28] += double(Count*Count);
    timer[27] += (endtime-starttime)*(endtime-starttime);
#endif /* MPI_INSTRUMENTATION */
  
  } // end: if (ProcessorNumber != ToProcessor)

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
    
    // Zero out NumberOfActiveParticles in the first pass, so we can
    // use grid::AddActiveParticles.
    if (type == 0) ToGrid->NumberOfActiveParticles = 0;
    
    // Extract the number of particles in the buffer
    // type_count should have been received by now...
    ap_info->UnpackBuffer(buffer, 0,
			  NewParticles, type_count[type]);

    for (i = 0; i < type_count[type]; i++)
    {
      NewParticles[i]->AssignCurrentGrid(ToGrid);
    }

    ToGrid->AddActiveParticles(NewParticles, 0, type_count[type]);
    
    delete[] buffer;
    
  } // end: if (MyProcessorNumber == ToProcessor && ...)
  
  } // ENDFOR particle types
  
  delete [] type_element_size;

#endif /* USE_MPI */
  return SUCCESS;
}

