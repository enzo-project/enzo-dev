/***********************************************************************
/
/  INLINE HALO FINDER HELPER FUNCTION :: FINALIZE
/
/  written by: John Wise
/  date:       June, 2009
/
/  PURPOSE:    Moves particle data from the halo finder's memory to the
/              enzo's memory.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
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

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

/************************************************************************/
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int CommunicationCollectParticles(LevelHierarchyEntry *LevelArray[],
				  int level, bool ParticlesAreLocal, 
				  int CollectMode);
/************************************************************************/

void FOF_Finalize(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
		  FOFData &D)
{
#ifdef UNUSED
  /* Get enzo units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, MassUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, MetaData->Time);
  
  // Sometimes MassUnits is infinite (in cgs) when using single
  // precision, so calculate it in double precision.

  double EnzoMassUnits = (double) DensityUnits * pow(LengthUnits, 3.0) /
    (MetaData->TopGridDims[0] * MetaData->TopGridDims[1] * 
     MetaData->TopGridDims[2]);

  /******************** MOVE PARTICLES BACK TO ENZO ********************/

  int i, level, lvl, gi, npart;
  LevelHierarchyEntry *Temp;
  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int nGrids[MAX_DEPTH_OF_HIERARCHY];
  
  /* Create grid lists */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    if (LevelArray[level] != NULL)
      nGrids[level] = GenerateGridArray(LevelArray, level, &Grids[level]);

  /* Count particles and then allocate memory for particles in grids */

  for (i = 1; i <= D.Nlocal; i++) {
    lvl = D.P[i].level;
    gi = D.P[i].GridID;
    npart = Grids[lvl][gi]->GridData->ReturnNumberOfParticles();
    Grids[lvl][gi]->GridData->SetNumberOfParticles(npart+1);
  }
  
  // Allocate memory, but reset NumberOfParticles to zero because
  // we'll use it as a counter when adding the particles back.

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
      npart = Temp->GridData->ReturnNumberOfParticles();
      if (npart > 0) {
	Temp->GridData->AllocateNewParticles(npart);
	Temp->GridData->SetNumberOfParticles(0);
      }
    } // ENDFOR grids

  /* Move particles back */

  for (i = 1; i <= D.Nlocal; i++) {
    lvl = D.P[i].level;
    gi = D.P[i].GridID;
    Grids[lvl][gi]->GridData->MoveParticlesFOF
      (lvl, gi, D.P, i, D, VelocityUnits, EnzoMassUnits, COPY_IN);
  } // ENDFOR particles

  /* Above we just copied the particles to the local grid, now we have
     to collect the particles on the grid's host processor
     (ProcessorNumber) */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    if (LevelArray[level] != NULL)
      CommunicationCollectParticles(LevelArray, level, false, SIBLINGS_ONLY);

  /* Clean up */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    if (LevelArray[level] != NULL)
      delete [] Grids[level];

#endif /* UNUSED */
  return;

}
