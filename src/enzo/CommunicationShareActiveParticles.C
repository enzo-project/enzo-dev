/***********************************************************************
/
/  COMMUNICATION ROUTINE: DISTRIBUTE ACTIVE PARTICLES TO PROCESSORS
/
/  written by: John Wise
/  date:       May, 2009
/  modified1   July, 2009 by John Wise: adapted for stars
/  modified2:  December, 2011 by John Wise: adapted for active particles
/  modified3:  November, 2012 by Nathan Goldbaum: Rewriting to use an Allgatherv
/
/  PURPOSE: Takes a list of active particles moves and sends/receives
/           them to all processors
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
#include "ActiveParticle.h"
#include "CommunicationUtilities.h"
#include "SortCompareFunctions.h"

void my_exit(int status);

int CommunicationShareActiveParticles(int *NumberToMove, 
				      ActiveParticleList<ActiveParticleType> &SendList,
				      int &NumberOfReceives, 
				      ActiveParticleList<ActiveParticleType> &SharedList)
{

  int i, type, proc, ap_id;
  int SizeOfSends, NumberOfNewParticles;
  ActiveParticleType_info *ap_info;

  int TotalNumberToMove = 0, GlobalNumberToMove;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];
  GlobalNumberToMove = TotalNumberToMove;
#ifdef USE_MPI
  CommunicationAllReduceValues(&GlobalNumberToMove, 1, MPI_SUM);
#endif
  if (GlobalNumberToMove == 0) return SUCCESS;

#ifdef USE_MPI
  MPI_Arg Count;
  MPI_Arg SendCount;
  MPI_Arg RecvCount;
  MPI_Arg stat;
#endif /* USE_MPI */

  if (NumberOfProcessors > 1) {

#ifdef USE_MPI

    /* We will have one collective communication per particle type.
       In the future, there might be a way to consolidate all of the
       buffers into one buffer and communicate it as a whole. */

    for (type = 0; type < EnabledActiveParticlesCount; type++) {

      ap_info = EnabledActiveParticles[type];
      ap_id = ap_info->GetEnabledParticleID();

      /* Get counts from each processor to allocate buffers. */

      int NumberToSend, local_buffer_size, instance_size;

      NumberToSend = 0;
      for (i = 0; i < TotalNumberToMove; i++)
        if (SendList[i]->ReturnType() == type) NumberToSend++;

      int *nCount = new int[NumberOfProcessors];

      MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

      MPI_Allgather(&NumberToSend, 1, DataTypeInt,
		    nCount, 1, DataTypeInt, MPI_COMM_WORLD);

      NumberOfNewParticles = 0;
      for (i = 0; i < NumberOfProcessors; i++) {
	NumberOfNewParticles += nCount[i];
      }

      /* Create a MPI packed buffer from the SendList */
      
      instance_size = ap_info->ReturnElementSize();

      if (NumberToSend != 0)
	local_buffer_size = NumberToSend*instance_size;
      else
	local_buffer_size = 0;

      int send_buffer_size;
      char *send_buffer = new char[local_buffer_size];
      // This will break if more than one AP type is in the simulation
      send_buffer_size = ap_info->FillBuffer(SendList, NumberToSend, send_buffer);

      /* generate displacement list. */

      Eint32 *displace = new Eint32[NumberOfProcessors];
      Eint32 *all_buffer_sizes = new Eint32[NumberOfProcessors];
      Eint32 total_buffer_size = 0, position = 0;

      for (i = 0; i < NumberOfProcessors; i++) {
	all_buffer_sizes[i] = nCount[i]*instance_size;
	total_buffer_size += all_buffer_sizes[i];
      }

      displace[0] = position;
      for (i = 1; i < NumberOfProcessors; i++) {
	if (all_buffer_sizes[i-1] > 0)
	  position += all_buffer_sizes[i-1];
	displace[i] = position;
      }

      position = 0;

      /* Allocate the receive buffer and gather the particle buffers */

      char* recv_buffer = new char[total_buffer_size];

      MPI_Allgatherv(send_buffer, local_buffer_size, MPI_PACKED,
		     recv_buffer, all_buffer_sizes, displace, MPI_PACKED, MPI_COMM_WORLD);

      /* Unpack the particle buffers, generate global shared active particle list */
      
      int count = 0;
      for (proc = 0; proc < NumberOfProcessors; proc++) {
        if (nCount[proc] > 0) {
          ap_info->UnpackBuffer(recv_buffer+displace[proc], count, 
              SharedList, nCount[proc]);
          count += nCount[proc];
        }
      }

      /* clean up */

      delete [] nCount;
      delete [] displace;
      delete [] send_buffer;
      delete [] recv_buffer;
      delete [] all_buffer_sizes;
      nCount = NULL;
      displace = NULL;
      send_buffer = NULL;
      recv_buffer = NULL;

      NumberOfReceives = NumberOfNewParticles;

#ifdef MPI_INSTRUMENTATION
      endtime = MPI_Wtime();
      timer[9] += endtime-starttime;
      counter[9] ++;
      timer[10] += double(NumberOfReceives);
      GlobalCommunication += endtime-starttime;
      CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
    } // ENDFOR types

#endif /* USE_MPI */    

  } // ENDIF multi-processor
  else {
    NumberOfReceives = TotalNumberToMove;
    SharedList = SendList;
  }

  // First sort the list by destination grid, so the searching for
  // grids is more efficient.

  SharedList.sort_grid(0, NumberOfReceives);

  return SUCCESS;

}
