/***********************************************************************
/
/  GRID CLASS (SORT PARTICLES BY PARTICLE TYPE)
/
/  written by: Greg Bryan
/  date:       Jan, 2001
/  modified1:  Matthew Turk
/  modified2:  John Wise (May 2010) -- Use shell sort
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
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
 
/* function prototypes */
 
void ShellSortAndDrag(int List[], int N,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, PINT  *DragList3[]);
 
 
void grid::SortParticlesByType()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfParticles == 0)
    return;
 
  int dim, j;
 
  /* Allocate arrays of pointer, one for float type and one for FLOAT type,
     and file them up with pointers to the particle data. */
 
  float **DragList1 = new float*[GridRank+1+NumberOfParticleAttributes];
  FLOAT **DragList2 = new FLOAT*[GridRank];
  PINT   **DragList3 = new PINT*[1];
  for (dim = 0; dim < GridRank; dim++) {
    DragList2[dim] = ParticlePosition[dim];
    DragList1[dim] = ParticleVelocity[dim];
  }
  DragList1[GridRank] = ParticleMass;
  DragList3[0]        = ParticleNumber;
  for (j = 0; j < NumberOfParticleAttributes; j++)
    DragList1[GridRank+1+j] = ParticleAttribute[j];
 
  /* Sort by particle index, dragging the data along. */
 
  ShellSortAndDrag(ParticleType, NumberOfParticles,
		   GridRank+1+NumberOfParticleAttributes, DragList1,
		   GridRank, DragList2, 1, DragList3);
 
  /* Clean up. */
 
  delete [] DragList1;
  delete [] DragList2;
  delete [] DragList3;
 
  return;
}
