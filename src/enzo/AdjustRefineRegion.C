/***********************************************************************
/
/  ADJUST THE REFINE REGION TO ONLY INCLUDE HIGH-RESOLUTION PARTICLES
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
#include "Star.h"
#include "CommunicationUtilities.h"

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData, int EL_level)
{

//if (!(RefineRegionAutoAdjust && EL_level == 0))
  if (!(RefineRegionAutoAdjust >= 1 && EL_level == RefineRegionAutoAdjust-1)) 
    return SUCCESS;

  if (RefineRegionLeftEdge[0] == DomainLeftEdge[0] &&
      RefineRegionRightEdge[0] == DomainRightEdge[0])
    return SUCCESS;

  int i, dim, idim, level;
  LevelHierarchyEntry *Temp;
  float MinParticleMass, dx;
  FLOAT RefineRegionWidth[MAX_DIMENSION], RefineRegionCenter[MAX_DIMENSION];
  int MaximumStaticRegionLevel = 0;

  /* First find the highest resolution particle in the whole
     simulation.  We could be clever by only looking in the finest
     levels, but there might be some false collapses caused by a
     massive particle entering the refine region.  */

  MinParticleMass = huge_number;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->FindMinimumParticleMass(MinParticleMass, level) == FAIL) {
	ENZO_FAIL("Error in grid::FindMinimumParticleMass.\n");
      }
  MinParticleMass = CommunicationMinValue(MinParticleMass);

  /* 
     Now we look for the bounding box that contains ONLY these high
     resolution particles.  To do this in a few passes through the
     data, it is more tricky than it seems at first glance.  It's not
     a simple min/max of the particle positions.  

     First we search for massive particles in the outer slab of cells
     on each face.  If there are any particles there, we adjust the
     region by one cell width in that dimension.

     It gets tricky when there are massive particles left that aren't in
     the outer slab but exist a few cell widths inside the region.
     This shouldn't happen often, but we have to account for it.  We
     look for this when no particles were removed in the last pass
     over all faces, but there are massive particles leftover.  In
     this case, we decrease the region width by one cell width until
     we've removed all of them.  

     To avoid favoring reducing the width on the first faces
     (x-dimension, in particular), I've randomized how the faces are
     looped.
  */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RefineRegionCenter[dim] = 0.5 * (RefineRegionLeftEdge[dim] + 
				     RefineRegionRightEdge[dim]);
    RefineRegionWidth[dim] = RefineRegionRightEdge[dim] - RefineRegionLeftEdge[dim];
  }

  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    MaximumStaticRegionLevel = max(MaximumStaticRegionLevel,
				   StaticRefineRegionLevel[i]+1);
  dx = (DomainRightEdge[0] - DomainLeftEdge[0]) / 
    (MetaData->TopGridDims[0] * pow(RefineBy, MaximumStaticRegionLevel));

  /* Make a local list of massive particles inside the refine region. */

  int face, side, dim1, dim2, NumberOfParticles = 0;
  int TotalNumberOfParticles, ParticlesLeft;
  FLOAT *ParticlePos[MAX_DIMENSION];

  // Count number of massive particles first, then allocate memory
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->FindMassiveParticles(MinParticleMass, level, 
				ParticlePos, NumberOfParticles, TRUE) == FAIL) {
	ENZO_FAIL("Error in grid::FindMassiveParticles(count).\n");
      }

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    ParticlePos[dim] = new FLOAT[NumberOfParticles];

  NumberOfParticles = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->FindMassiveParticles(MinParticleMass, level, 
				ParticlePos, NumberOfParticles, FALSE) == FAIL) {
	ENZO_FAIL("Error in grid::FindMassiveParticles.\n");
      }

  // Define some convenient variables, such as (1) a flagging field
  // that indicates which particles have been removed, and (2) an
  // array with particle positions in units of cell widths.

  bool inside;
  bool* RemoveFlag = new bool[NumberOfParticles];
  int* CellPosition[MAX_DIMENSION];
  int RefineRegionLeftEdgeCell[MAX_DIMENSION], RefineRegionRightEdgeCell[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    CellPosition[dim] = new int[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      CellPosition[dim][i] = (int) (ParticlePos[dim][i] / dx);
    delete [] ParticlePos[dim];
  }

  for (i = 0; i < NumberOfParticles; i++)
    RemoveFlag[i] = false;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RefineRegionLeftEdgeCell[dim]  = (int) (RefineRegionLeftEdge[dim] / dx);
    RefineRegionRightEdgeCell[dim] = (int) (RefineRegionRightEdge[dim] / dx);
  }
  
  // Refine boundaries until no massive particles are included inside.
  // Do this in units of cell widths since it's easier (and faster) to
  // deal with integers.
  
  TotalNumberOfParticles = NumberOfParticles;
#ifdef USE_MPI
  CommunicationAllReduceValues(&TotalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */

  int irand, threshold;
  int nRemoveLast, nRemoveTotal, nRemoveCurrent, NumberToRemove;

  nRemoveTotal = INT_UNDEFINED;
  ParticlesLeft = TotalNumberOfParticles;
  if (rand_init == 0) {
    srand( time(NULL) );
    rand_init = 1;
  }

  while (ParticlesLeft > 0) {

    nRemoveLast = nRemoveTotal;
    nRemoveTotal = 0;

    for (irand = 0; irand < 6; irand++) {

      // Randomize which face we pick in order to not favor refining
      // the first (x-dim) ones
      if (MyProcessorNumber == ROOT_PROCESSOR)
	face = rand() % 6;
      CommunicationBroadcastValue(&face, ROOT_PROCESSOR);
      dim = face/2;
      side = face%2;
      dim1 = (dim+1) % MAX_DIMENSION;
      dim2 = (dim+2) % MAX_DIMENSION;
      
      // Look for particles in the first and last slab of cells,
      // excluding the cells on edges.

      NumberToRemove = 0;
      if (side == 0) {

	// Left face

	threshold = RefineRegionLeftEdgeCell[dim] + 1;
	for (i = 0; i < NumberOfParticles; i++)
	  if (!RemoveFlag[i])
	    if (CellPosition[dim][i] <= threshold &&
		CellPosition[dim1][i] > RefineRegionLeftEdgeCell[dim1] &&
		CellPosition[dim1][i] < RefineRegionRightEdgeCell[dim1] &&
		CellPosition[dim2][i] > RefineRegionLeftEdgeCell[dim2] &&
		CellPosition[dim2][i] < RefineRegionRightEdgeCell[dim2])
	      NumberToRemove++;

	if (NumberToRemove > 0 || nRemoveLast == 0) {
	  nRemoveCurrent = 0;
	  RefineRegionLeftEdgeCell[dim]++;
	  for (i = 0; i < NumberOfParticles; i++)
	    if (!RemoveFlag[i]) {
	      inside = true;
	      for (idim = 0; idim < MAX_DIMENSION; idim++)
		inside &= (CellPosition[idim][i] > RefineRegionLeftEdgeCell[idim] &&
			   CellPosition[idim][i] < RefineRegionRightEdgeCell[idim]);
	      if (!inside) {
		RemoveFlag[i] = true;
		nRemoveCurrent++;
	      }
	    } // ENDIF !RemoveFlag

	  // If we've removed any particles in the case where we
	  // didn't remove any particles in the previous pass
	  // (nRemoveLast == 0), reset the counter.
	  if (nRemoveCurrent > 0)
	    nRemoveLast = INT_UNDEFINED;

	} // ENDIF any particles to remove

      } // ENDIF side == 0
      else {

	// Right face

	threshold = RefineRegionRightEdgeCell[dim] - 1;
	for (i = 0; i < NumberOfParticles; i++)
	  if (!RemoveFlag[i])
	    if (CellPosition[dim][i] >= threshold &&
		CellPosition[dim1][i] > RefineRegionLeftEdgeCell[dim1] &&
		CellPosition[dim1][i] < RefineRegionRightEdgeCell[dim1] &&
		CellPosition[dim2][i] > RefineRegionLeftEdgeCell[dim2] &&
		CellPosition[dim2][i] < RefineRegionRightEdgeCell[dim2])
	      NumberToRemove++;

	if (NumberToRemove > 0 || nRemoveLast == 0) {
	  nRemoveCurrent = 0;
	  RefineRegionRightEdgeCell[dim]--;
	  for (i = 0; i < NumberOfParticles; i++)
	    if (!RemoveFlag[i]) {
	      inside = true;
	      for (idim = 0; idim < MAX_DIMENSION; idim++)
		inside &= (CellPosition[idim][i] > RefineRegionLeftEdgeCell[idim] &&
			   CellPosition[idim][i] < RefineRegionRightEdgeCell[idim]);
	      if (!inside) {
		RemoveFlag[i] = true;
		nRemoveCurrent++;
	      }
	    } // ENDIF !RemoveFlag

	  // If we've removed any particles in the case where we
	  // didn't remove any particles in the previous pass
	  // (nRemoveLast == 0), reset the counter.
	  if (nRemoveCurrent > 0)
	    nRemoveLast = INT_UNDEFINED;

	} // ENDIF any particles to remove

      } // ENDELSE side == 0

      nRemoveTotal += NumberToRemove;
      ParticlesLeft = 0;
      for (i = 0; i < NumberOfParticles; i++)
	if (!RemoveFlag[i]) ParticlesLeft++;

      // Synchronize over all processors
#ifdef USE_MPI
      CommunicationAllReduceValues(&ParticlesLeft, 1, MPI_SUM);
      CommunicationAllReduceValues(&nRemoveLast, 1, MPI_MIN);
      CommunicationAllReduceValues(RefineRegionLeftEdgeCell, 
				   MAX_DIMENSION, MPI_MAX);
      CommunicationAllReduceValues(RefineRegionRightEdgeCell, 
				   MAX_DIMENSION, MPI_MIN);
#endif /* USE_MPI */

//      if (debug) {
//	printf("LeftEdge = %"ISYM" %"ISYM" %"ISYM"\n", RefineRegionLeftEdgeCell[0], 
//	       RefineRegionLeftEdgeCell[1], RefineRegionLeftEdgeCell[2]);
//	printf("RightEdge = %"ISYM" %"ISYM" %"ISYM"\n", RefineRegionRightEdgeCell[0], 
//	       RefineRegionRightEdgeCell[1], RefineRegionRightEdgeCell[2]);
//	printf("dim %"ISYM", side %"ISYM" :: removed %"ISYM", %"ISYM" particles left\n",
//	       dim, side, NumberToRemove, ParticlesLeft);
//      }

      if (ParticlesLeft == 0)
	break;

      for (dim = 0; dim < MAX_DIMENSION; dim++)
	if (RefineRegionLeftEdgeCell[dim] >= RefineRegionRightEdgeCell[dim]) {
	  fprintf(stderr, "Refine region collapsed to nothing!\n");
	  fprintf(stderr, "RefineRegionLeftEdgeCell = %"ISYM" %"ISYM" %"ISYM"\n", 
		  RefineRegionLeftEdgeCell[0], RefineRegionLeftEdgeCell[1], 
		  RefineRegionLeftEdgeCell[2]);
	  fprintf(stderr, "RefineRegionRightEdgeCell = %"ISYM" %"ISYM" %"ISYM"\n", 
		  RefineRegionRightEdgeCell[0], RefineRegionRightEdgeCell[1], 
		  RefineRegionRightEdgeCell[2]);
	  ENZO_FAIL("Refine region collapsed to nothing!\n");
	}

    } // ENDFOR region faces

    // Synchronize over all processors
#ifdef USE_MPI
    CommunicationAllReduceValues(&nRemoveTotal, 1, MPI_SUM);
#endif /* USE_MPI */

  } // ENDWHILE particles left

  // Adjust the RefineRegion with the new cell-based region
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RefineRegionLeftEdge[dim] = RefineRegionLeftEdgeCell[dim] * dx;
    RefineRegionRightEdge[dim] = RefineRegionRightEdgeCell[dim] * dx;
  }
  
  if (MyProcessorNumber == ROOT_PROCESSOR && TotalNumberOfParticles > 0) {

    printf("AdjustRefineRegion: Changed RefineRegionLeftEdge to "
	   "[%"FSYM" %"FSYM" %"FSYM"]\n", 
	   RefineRegionLeftEdge[0], RefineRegionLeftEdge[1], 
	   RefineRegionLeftEdge[2]);
    printf("AdjustRefineRegion: Changed RefineRegionRightEdge to "
	   "[%"FSYM" %"FSYM" %"FSYM"]\n", 
	   RefineRegionRightEdge[0], RefineRegionRightEdge[1], 
	   RefineRegionRightEdge[2]);
  }

  /* Clean up */

  delete [] RemoveFlag;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    delete [] CellPosition[dim];

  return SUCCESS;

}
