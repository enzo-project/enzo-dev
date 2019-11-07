/***********************************************************************
/
/  COMMUNICATION ROUTINE: SYNCHRONIZE GLOBAL ACTIVE PARTICLE LIST ACROSS
/                         PROCESSORS
/
/  written by: John Wise
/  date:       March 2009
/  modified1:  Nathan Goldbaum, March 2012
/              Porting to active particles
/
/
/  PURPOSE: Generates a list of all active particles in the simulation
/           with id = ActiveParticleIDToFind.
/
*************************************************************************/

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
#include "SortCompareFunctions.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
void ActiveParticleFindAll(
    LevelHierarchyEntry *LevelArray[], int *GlobalNumberOfActiveParticles, 
    int ActiveParticleIDToFind, 
    ActiveParticleList<ActiveParticleType>& GlobalList)
{
  int i, level, type, ap_id, GridNum, LocalNumberOfActiveParticles, proc, 
    buffer_size, LocalNumberOfActiveParticlesOnThisLevel, element_size, count, 
    offset;
  ActiveParticleList<ActiveParticleType> LocalActiveParticlesOfThisType, 
    LocalActiveParticlesOnThisLevel, ParticlesOnThisProc;

  HierarchyEntry **Grids = NULL;
  int NumberOfGrids, *NumberOfActiveParticlesInGrids = NULL;
  ActiveParticleType_info *ap_info = NULL;

  *GlobalNumberOfActiveParticles = 0;
  LocalNumberOfActiveParticles = 0;

  for (type = 0; type < EnabledActiveParticlesCount; type++) {
      
    ap_info = EnabledActiveParticles[type];
    ap_id = ap_info->GetEnabledParticleID();
    
    if (ap_id == ActiveParticleIDToFind) {
      
      /* Traverse the hierarchy and generate a buffer of all of the active 
	 particles of this type on this processor */
      
      for (level = 0; level <= MaximumRefinementLevel; level++) {
	NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
	NumberOfActiveParticlesInGrids = new int[NumberOfGrids];
	LocalNumberOfActiveParticlesOnThisLevel = 0;

	/* In a first pass, find the number of active particles on each grid */
	for(GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
	  
	  NumberOfActiveParticlesInGrids[GridNum] = Grids[GridNum]->GridData->
	    ReturnNumberOfActiveParticlesOfThisType(ActiveParticleIDToFind);
	  LocalNumberOfActiveParticlesOnThisLevel += NumberOfActiveParticlesInGrids[GridNum];
	} /* ENDFOR grids */
	
	offset = 0;
	
	if (LocalNumberOfActiveParticlesOnThisLevel > 0) {

	  /* If we've already found active particles, save the list */

	  if (LocalNumberOfActiveParticles != 0) 
	    {
	      ParticlesOnThisProc = LocalActiveParticlesOfThisType;
	      LocalActiveParticlesOfThisType.clear();
	    }
      
      /* In a second pass, fill up the active particle list for this level*/
	  for(GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
        // I've removed this function from the grid class, do this another way.
	    Grids[GridNum]->GridData->
	      AppendActiveParticlesToList(LocalActiveParticlesOnThisLevel,
              ActiveParticleIDToFind);
	    offset += NumberOfActiveParticlesInGrids[GridNum];
	  } 

	  LocalActiveParticlesOfThisType.reserve(
          offset + LocalNumberOfActiveParticles);
	  
	  /* If we've already found active particles, copy the cached
	     list to the new one and delete the old list */
	  if (LocalNumberOfActiveParticles != 0) {
	    for (i = 0; i < LocalNumberOfActiveParticles; i++)
	      LocalActiveParticlesOfThisType.copy_and_insert(
              *ParticlesOnThisProc[i]);
	    ParticlesOnThisProc.clear();
	  }

	  /* Finally, append the new active particles to the list */
	  for(i = LocalNumberOfActiveParticles; i < offset + LocalNumberOfActiveParticles; i++) 
	    LocalActiveParticlesOfThisType.copy_and_insert(
            *LocalActiveParticlesOnThisLevel[i - LocalNumberOfActiveParticles]);

	  LocalNumberOfActiveParticles += LocalNumberOfActiveParticlesOnThisLevel;
	    
	  /* Delete the list for this level */
      LocalActiveParticlesOnThisLevel.clear();

	} 
	
	delete [] Grids;
	Grids = NULL;
	delete [] NumberOfActiveParticlesInGrids;
	NumberOfActiveParticlesInGrids = NULL;

      } /* ENDFOR level */
      
    } /* ENDIF id == id to search for */

    /**************************************************/
    /*                                                */
    /* Share active particle counts on all processors */
    /*                                                */
    /**************************************************/

    int *nCount = NULL;
        
    if (NumberOfProcessors > 1) {
#ifdef USE_MPI
      
      nCount = new int[NumberOfProcessors];
      
      MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

      MPI_Allgather(&LocalNumberOfActiveParticles, 1, DataTypeInt, 
		    nCount, 1, DataTypeInt, MPI_COMM_WORLD);
      
      for (i = 0; i < NumberOfProcessors; i++) {
	*GlobalNumberOfActiveParticles += nCount[i];
      }
      
#endif /* USE_MPI */
    } /* ENDIF Number of processors > 1 */
    else {
      *GlobalNumberOfActiveParticles = LocalNumberOfActiveParticles;
    }
    /**************************************************/
    /*                                                */
    /* Gather the active particles on all processors  */
    /*                                                */
    /**************************************************/
    
    if (*GlobalNumberOfActiveParticles > 0) {

      Eint32 *displace = new Eint32[NumberOfProcessors];
      Eint32 *all_buffer_sizes = new Eint32[NumberOfProcessors];

      GlobalList.reserve(*GlobalNumberOfActiveParticles);

      if (NumberOfProcessors > 1) {
	
#ifdef USE_MPI
	/* Construct the MPI packed  buffer from the list of local particles*/
	Eint32 total_buffer_size=0, local_buffer_size, position = 0;
	char *send_buffer = NULL, *recv_buffer = NULL;
	element_size = ap_info->ReturnElementSize();
	
	for (i = 0; i < NumberOfProcessors; i++) {
	  all_buffer_sizes[i] = nCount[i]*element_size;
	  total_buffer_size += all_buffer_sizes[i];
	}

	displace[0] = position;
	for (i = 1; i < NumberOfProcessors; i++) {
	  if (nCount[i-1] > 0)
	    position += all_buffer_sizes[i-1];
	  displace[i] = position;
	}

	position = 0;

	if (LocalNumberOfActiveParticles > 0)
	  local_buffer_size = LocalNumberOfActiveParticles*element_size;
	else
	  local_buffer_size = 0;

	send_buffer = new char[local_buffer_size];
	recv_buffer = new char[total_buffer_size];
	
	ap_info->FillBuffer(LocalActiveParticlesOfThisType,
                        LocalNumberOfActiveParticles,
                        send_buffer);

	/* Share all data with all processors */

	MPI_Allgatherv(send_buffer, local_buffer_size, MPI_PACKED,
		       recv_buffer, all_buffer_sizes, displace, MPI_PACKED, MPI_COMM_WORLD);

	/* Unpack MPI buffers, generate global active particles list */

	count = 0;
	for (proc = 0; proc < NumberOfProcessors; proc++) {
	  if (nCount[proc] > 0) {
	    ap_info->UnpackBuffer(recv_buffer+displace[proc], count,
				  GlobalList, nCount[proc]);
	    count += nCount[proc];
	  }
	}

	delete [] nCount;
	delete [] displace;
	delete [] send_buffer;
	delete [] recv_buffer;
	delete [] all_buffer_sizes;
	nCount = NULL;
	displace = NULL;
	send_buffer = NULL;
	recv_buffer = NULL;

	/* Set grid pointers */
	for (i = 0; i < *GlobalNumberOfActiveParticles; i++) {
	  NumberOfGrids = GenerateGridArray(LevelArray, 
          GlobalList[i]->ReturnLevel(), &Grids);
	  for (GridNum = 0; GridNum < NumberOfGrids; GridNum++) {
	    if (Grids[GridNum]->GridData->GetGridID() == 
            GlobalList[i]->ReturnGridID())
	      GlobalList[i]->AssignCurrentGrid(Grids[GridNum]->GridData);
	  }
	  delete [] Grids;
	}

	LocalActiveParticlesOfThisType.clear();

#endif /* USE_MPI */
       
      } /* ENDIF multi-processor */
      else {
        GlobalList = LocalActiveParticlesOfThisType;
        LocalActiveParticlesOfThisType.clear();
      } // ENDIF serial
      
    }  /* ENDIF number of active particles > 0 */
    else 
      delete [] nCount;
    
  } /* ENFOR Active particle types */

}
  
