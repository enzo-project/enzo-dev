/***********************************************************************
/
/  GRID CLASS (COPY SUBGRID PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: modified to move subgrid particles
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
 
int grid::TransferSubgridParticles(grid* Subgrids[], int NumberOfSubgrids, 
				   int* &NumberToMove, int StartIndex, 
				   int EndIndex, particle_data* &List, 
				   bool KeepLocal, bool ParticlesAreLocal,
				   int CopyDirection, int IncludeGhostZones,
				   int CountOnly)
{
 
  /* Declarations. */

  int i, j, index, dim, n1, grid, proc;
  int i0, j0, k0;

  /* ----------------------------------------------------------------- */
  /* Copy particle out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If particles aren't distributed over several processors, exit
       if this isn't the host processor. */

    if (ParticlesAreLocal && MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* If there are no particles to move, we're done. */

    if (NumberOfParticles == 0) {
      delete [] BaryonField[NumberOfBaryonFields];
      BaryonField[NumberOfBaryonFields] = NULL;
      return SUCCESS;
    }

    /* Count the number of particles already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];

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
 
    /* Count particles to move */

    int *subgrid = NULL;
    subgrid = new int[NumberOfParticles];

    for (i = 0; i < NumberOfParticles; i++) {

      /* Compute index of particle position. */
 
      i0 = int((ParticlePosition[0][i] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((ParticlePosition[1][i] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((ParticlePosition[2][i] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(EndIndex[0], i0), StartIndex[0]);
      j0 = max(min(EndIndex[1], j0), StartIndex[1]);
      k0 = max(min(EndIndex[2], k0), StartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find and store subgrid number of this particle, and add to
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
	ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", 
		subgrid[i], NumberOfSubgrids)
      }
      
    } // ENDFOR particles

    if (CountOnly == TRUE) {
      delete [] subgrid;
      return SUCCESS;
    }

    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];

    if (TotalToMove > PreviousTotalToMove) {

      // Increase the size of the list to include the particles from
      // this grid

//      particle_data *NewList = new particle_data[TotalToMove];
//      memcpy(NewList, List, PreviousTotalToMove * sizeof(particle_data));
//      delete [] List;
//      List = NewList;
//      particle_data *TempList = List;
//      List = new particle_data[TotalToMove];
//      for (i = 0; i < PreviousTotalToMove; i++)
//	List[i] = TempList[i];
//      delete [] TempList;

      /* Compute the increase in mass for particles moving to the subgrid. */
 
      float RefinementFactors[MAX_DIMENSION], MassIncrease = 1.0;
      this->ComputeRefinementFactorsFloat(Subgrids[0], RefinementFactors);
      for (dim = 0; dim < GridRank; dim++)
	MassIncrease *= RefinementFactors[dim];

      /* Move particles (mark those moved by setting mass = FLOAT_UNDEFINED). */

      n1 = PreviousTotalToMove;

      for (i = 0; i < NumberOfParticles; i++) {
	if (subgrid[i] >= 0) {
	  if (KeepLocal)
	    proc = MyProcessorNumber;
	  else
	    proc = Subgrids[subgrid[i]]->ReturnProcessorNumber();
	  for (dim = 0; dim < GridRank; dim++) {
	    List[n1].pos[dim] = ParticlePosition[dim][i];
	    List[n1].vel[dim] = ParticleVelocity[dim][i];
	  }
	  List[n1].mass = ParticleMass[i] * MassIncrease;
	  List[n1].id = ParticleNumber[i];
	  List[n1].type = ParticleType[i];
	  for (j = 0; j < NumberOfParticleAttributes; j++)
	    List[n1].attribute[j] = ParticleAttribute[j][i];
	  List[n1].grid = subgrid[i];
	  List[n1].proc = proc;
	  ParticleMass[i] = FLOAT_UNDEFINED;
	  n1++;
	} // ENDIF move to subgrid
      } // ENDFOR particles

      if (TotalToMove != PreviousTotalToMove)
	this->CleanUpMovedParticles();

    } // ENDIF particles to move

    delete [] subgrid;
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy particle back into grid. */
 
  else {

    /* If particles aren't distributed over several processors, exit
       if this isn't the host processor. */

    if (ParticlesAreLocal && MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */
 
    int TotalNumberOfParticles;
    int NumberOfNewParticles = EndIndex - StartIndex;

    TotalNumberOfParticles = NumberOfParticles + NumberOfNewParticles;

    /* Allocate space for the particles. */

    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION], *Mass,
          *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    PINT *Number;
    int  *Type;
 
    if (TotalNumberOfParticles > 0) {

    Mass = new float[TotalNumberOfParticles];
    Number = new PINT[TotalNumberOfParticles];
    Type = new int[TotalNumberOfParticles];
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim] = new FLOAT[TotalNumberOfParticles];
      Velocity[dim] = new float[TotalNumberOfParticles];
    }
    for (i = 0; i < NumberOfParticleAttributes; i++)
      Attribute[i] = new float[TotalNumberOfParticles];

    if (Velocity[GridRank-1] == NULL && TotalNumberOfParticles != 0) {
      fprintf(stderr, "malloc error (out of memory?)\n");
      ENZO_FAIL("Error in Grid_TransferSubgridParticles.C");
    }

    /* Copy this grid's particles to the new space. */

    for (i = 0; i < NumberOfParticles; i++) {
      Mass[i] = ParticleMass[i];
      Number[i] = ParticleNumber[i];
      Type[i] = ParticleType[i];
    }

    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < NumberOfParticles; i++) {
	Position[dim][i] = ParticlePosition[dim][i];
	Velocity[dim][i] = ParticleVelocity[dim][i];
      }
	
    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < NumberOfParticles; i++)
	  Attribute[j][i] = ParticleAttribute[j][i];
 
    /* Copy new particles */

    int n = NumberOfParticles;

    for (i = StartIndex; i < EndIndex; i++) {
      Mass[n] = List[i].mass;
      Number[n] = List[i].id;
      Type[n] = List[i].type;
      n++;
    }

    for (dim = 0; dim < GridRank; dim++) {
      n = NumberOfParticles;
      for (i = StartIndex; i < EndIndex; i++) {
	Position[dim][n] = List[i].pos[dim];
	Velocity[dim][n] = List[i].vel[dim];
	n++;
      }
    }
      
    for (j = 0; j < NumberOfParticleAttributes; j++) {
      n = NumberOfParticles;
      for (i = StartIndex; i < EndIndex; i++) {
	Attribute[j][n] = List[i].attribute[j];
	n++;
      }
    }
      
    } // ENDIF TotalNumberOfParticles > 0

    /* Set new number of particles in this grid. */
 
    NumberOfParticles = TotalNumberOfParticles;
 
    /* Copy new ones into place. */
 
    this->DeleteParticles();
    if (TotalNumberOfParticles > 0)
      this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
				Attribute);
 
  } // end: if (COPY_IN)


 
  return SUCCESS;
}
