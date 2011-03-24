/***********************************************************************
/
/  INLINE HALO FINDER HELPER FUNCTION :: INITIALIZATION
/
/  written by: John Wise
/  date:       November, 2009
/
/  PURPOSE:    Moves particle data from FOF memory to the enzo's 
/              grid memory.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids,
				   int TopGridDims[]);
int CommunicationCollectParticles(LevelHierarchyEntry *LevelArray[],
				  int level, bool ParticlesAreLocal, 
				  bool SyncNumberOfParticles, 
				  bool MoveStars, int CollectMode);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);

/* 
   Move all particles on this processor (slab) to the first root grid
   on this processor, then let CommunicationTransferParticles and
   CommunicationCollectParticles do all the work.
*/

void FOF_Finalize(FOFData &D, LevelHierarchyEntry *LevelArray[], 
		  TopGridData *MetaData, int FOFOnly)
{

  /* Get enzo units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits;
  double MassUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, MetaData->Time);

  // Sometimes MassUnits is infinite (in cgs) when using single
  // precision, so calculate it in double precision.

  int TopGridDims3 = MetaData->TopGridDims[0] * MetaData->TopGridDims[1] * 
    MetaData->TopGridDims[2];

  MassUnits = (double) DensityUnits * pow(LengthUnits, 3.0) /
    TopGridDims3;

  int dummy, Zero = 0;
  LevelHierarchyEntry *Temp;

  D.P++;  // convert to zero-based for convenience
  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
    if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      Temp->GridData->MoveParticlesFOF(Zero, D.P, dummy, D, 
				       VelocityUnits, MassUnits, COPY_IN);
      break;
    }

  /* Particles are copied into the grid, now delete them in the FOF
     variable. */

  delete [] D.P;

  if (FOFOnly == FALSE) {

  /* Move the particles into their correct grids.  Just like
     RebuildHierarchy. */

  int level, n, ngrids = 0;
  grid *GridPointer[MAX_NUMBER_OF_TASKS];
  HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_TASKS];

  for (Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel)
    ngrids++;

  for (n = 0, Temp = LevelArray[0]; Temp; Temp = Temp->NextGridThisLevel) {
    GridPointer[n] = Temp->GridData;
    GridHierarchyPointer[n] = Temp->GridHierarchyEntry;
    n++;
  }

  bool ParticlesAreLocal, SyncNumberOfParticles, MoveStars;

  // Root grids (sync number of particles first)
  ParticlesAreLocal = false;
  SyncNumberOfParticles = false;
  MoveStars = false;
  CommunicationCollectParticles(LevelArray, 0, ParticlesAreLocal,
				SyncNumberOfParticles, MoveStars,
				SIBLINGS_ONLY);

  CommunicationTransferParticles(GridPointer, ngrids, MetaData->TopGridDims);

  ParticlesAreLocal = false;
  SyncNumberOfParticles = true;
  CommunicationCollectParticles(LevelArray, 0, ParticlesAreLocal, 
				SyncNumberOfParticles, MoveStars,
				SIBLINGS_ONLY);

  // Move to subgrids.  SUBGRIDS_GLOBAL because subgrids could be
  // distributed across many processors.
  ParticlesAreLocal = true;
  SyncNumberOfParticles = true;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    if (LevelArray[level] != NULL)
      CommunicationCollectParticles(LevelArray, level, ParticlesAreLocal,
				    SyncNumberOfParticles, MoveStars,
				    SUBGRIDS_GLOBAL);

  } // ENDIF FOFOnly == FALSE

}
