/***********************************************************************
/
/  GRID CLASS (COPY SUBGRID ACTIVE PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: modified to move subgrid particles
/  modified3:  July, 2009 by John Wise: modified to move stars
/  modified4:  December, 2011 by John Wise: modified for active particles
/
/  PURPOSE:
/
************************************************************************/
 
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"

int grid::TransferSubgridActiveParticles
(grid* Subgrids[], int NumberOfSubgrids, int* &NumberToMove, int StartIndex, 
 int EndIndex, ActiveParticleList<ActiveParticleType> &List, bool KeepLocal, 
 bool ParticlesAreLocal, int CopyDirection, int IncludeGhostZones, 
 int CountOnly)
{
 
  /* Declarations. */

  int i, j, index, dim, n1, grid, proc, type;
  int i0, j0, k0;

  /* ----------------------------------------------------------------- */
  /* Copy stars out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If stars aren't distributed over several processors, exit
       if this isn't the host processor. */

    if (ParticlesAreLocal && MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* If there are no active particles to move, we're done. */

    if (NumberOfActiveParticles == 0)
      return SUCCESS;

    /* Set boundaries (with and without ghost zones) */

    int StartIndex[] = {1,1,1}, EndIndex[] = {1,1,1};
    if (IncludeGhostZones)
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = 0;
        EndIndex[dim] = GridDimension[dim]-1;
      }
    else
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = GridStartIndex[dim];
        EndIndex[dim] = GridEndIndex[dim];
      }
 
    /* Count the number of stars already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];
 
    /* Count stars to move */

    int *subgrid = NULL;
    subgrid = new int[NumberOfActiveParticles];

    for (i = 0; i < NumberOfActiveParticles; i++) {

      /* Compute index of star position. */
 
      i0 = int((this->ActiveParticles[i]->pos[0] - CellLeftEdge[0][0]) / 
          CellWidth[0][0]);
      if (GridRank > 0)
        j0 = int((this->ActiveParticles[i]->pos[1] - CellLeftEdge[1][0]) / 
            CellWidth[1][0]);
      if (GridRank > 1)
        k0 = int((this->ActiveParticles[i]->pos[2] - CellLeftEdge[2][0]) / 
            CellWidth[2][0]);
 
      i0 = max(min(EndIndex[0], i0), StartIndex[0]);
      j0 = max(min(EndIndex[1], j0), StartIndex[1]);
      k0 = max(min(EndIndex[2], k0), StartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find and store subgrid number of this star, and add to
	 count. */

      subgrid[i] = nint(BaryonField[NumberOfBaryonFields][index])-1;
      if (subgrid[i] >= 0) {
        if (KeepLocal)
          proc = MyProcessorNumber;
        else
          proc = Subgrids[subgrid[i]]->ReturnProcessorNumber();
        NumberToMove[proc]++;
      }
      if (subgrid[i] < -1 || subgrid[i] > NumberOfSubgrids-1) {
        ENZO_VFAIL("star subgrid (%"ISYM"/%"ISYM") out of range\n", 
            subgrid[i], NumberOfSubgrids)
          }
      
    } // ENDFOR particles
    
    if (CountOnly == TRUE) {
      delete[] subgrid;
      return SUCCESS;
    }

    /* Allocate space. */

    int ParticlesLeft, NumberToMoveThisGrid;
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];

    if (TotalToMove > PreviousTotalToMove) {

      /* Compute the increase in mass for particles moving to the subgrid. */
 
      float RefinementFactors[MAX_DIMENSION], MassIncrease = 1.0;
      this->ComputeRefinementFactorsFloat(Subgrids[0], RefinementFactors);
      for (dim = 0; dim < GridRank; dim++)
        MassIncrease *= RefinementFactors[dim];

      /* Move active particles */
      
      NumberToMoveThisGrid = TotalToMove - PreviousTotalToMove;
      ParticlesLeft = NumberOfActiveParticles - NumberToMoveThisGrid;

      /* Move particles from grid array to a separate list. */

      n1 = PreviousTotalToMove;
      
      for (i = 0; i < NumberOfActiveParticles; i++) {
        if (subgrid[i] >= 0) {
          List.copy_and_insert(*ActiveParticles[i]);

          ActiveParticles.mark_for_deletion(i);

          List[n1]->SetGridID(subgrid[i]);
          // Increase the level if moving to a subgrid
          if (IncludeGhostZones == FALSE) {
            List[n1]->IncreaseLevel();
            List[n1]->AdjustMassByFactor(MassIncrease);
          }
          n1++;
        } // ENDIF subgrid
      } // ENDFOR particles

      ActiveParticles.delete_marked_particles();

      NumberOfActiveParticles = ParticlesLeft;

    } // ENDIF stars to move

    delete [] subgrid;
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    /* Count up total number. */
 
    int NumberOfNewActiveParticles = EndIndex - StartIndex;

    /* Copy stars from buffer into linked list */
    
    if (NumberOfNewActiveParticles > 0) {

      // Increase the level if moving to a subgrid
//      if (IncludeGhostZones == FALSE)
//	for (i = StartIndex; i < EndIndex; i++) {
//	}
      
      this->AddActiveParticles(List, StartIndex, EndIndex);

    } // ENDIF new particles

  } // end: if (COPY_IN)
 
  return SUCCESS;
}
