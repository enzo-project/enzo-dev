/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (HANDLE PARTICLE EXTERNAL BOUNDARY CONDITIONS)
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
//
//
int ExternalBoundary::SetExternalBoundaryParticles(int FieldRank,
						  int NumberOfParticles,
						  FLOAT *Position[],
						  float *Velocity[])
{
 
  /* declarations */
 
  int dim, i;
 
  /* Error check: grid ranks. */
 
  if (FieldRank != BoundaryRank) {
    ENZO_VFAIL("FieldRank(%"ISYM") != BoundaryRank(%"ISYM").\n",
            FieldRank, BoundaryRank)
  }
 
  /* Error check: allowed boundary types. */
 
  if (ParticleBoundaryType != periodic) {
    ENZO_FAIL("only periodic particle boundary conditions supported.\n");
  }
 
  /* PERIODIC BOUNDARY: wrap particles in each dimension. */
 
  if (ParticleBoundaryType == periodic && NumberOfProcessors == 1)
 
    for (dim = 0; dim < FieldRank; dim++)
 
      for (i = 0; i < NumberOfParticles; i++) {
 
	if (Position[dim][i] < DomainLeftEdge[dim])
	  Position[dim][i] += DomainRightEdge[dim] - DomainLeftEdge[dim];
 
	if (Position[dim][i] > DomainRightEdge[dim])

	  Position[dim][i] -= DomainRightEdge[dim] - DomainLeftEdge[dim];
 
      }	
 
  return SUCCESS;
 
}
