/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Alexei Kritsuk
/  date:       August, 2007
/
/  PURPOSE: Generates RandomForcingField[] arrays.
/           mode = 0 ==> fill out velocity fields
/           mode = 1 ==> fill out random force
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
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

extern "C" void FORTRAN_NAME(turboinit)(int *rank, int *nbox,
                          float *u, float *v, float *w,
                          int *in, int *jn, int *kn,
					int *ig, int *jg, int *kg, FLOAT * RandomMachNumber);

int grid::ComputeRandomForcingFields(int mode)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields < 3)
    ERROR_MESSAGE;


  /* initialize */

  int dim, i, size, nbox;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int Offset[MAX_DIMENSION];

  /* Compute size (in floats) of the current grid. */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* Find fields: density, total energy, velocity1-3. */
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Velocity1 must exist, but if 2 & 3 aren't present, then create blank
     buffers for them. */

  float *velocity1, *velocity2, *velocity3;
  if (mode == 0)
    velocity1 = BaryonField[Vel1Num];
  if (mode == 1) {
    if (RandomForcingField[0] == NULL)
      RandomForcingField[0] = new float[size];
    velocity1 = RandomForcingField[0];
  }
  if (mode < 0 || mode > 1)
    ERROR_MESSAGE;

  if (GridRank > 1)
    if (mode == 0)
      velocity2 = BaryonField[Vel2Num];
    else {
      if (RandomForcingField[1] == NULL)
	RandomForcingField[1] = new float[size];
      velocity2 = RandomForcingField[1];
    }
  else {
    velocity2 = new float[size];
    for (i = 0; i < size; i++)
      velocity2[i] = 0.0;
  }

  if (GridRank > 2)
    if (mode == 0)
      velocity3 = BaryonField[Vel3Num];
    else {
      if (RandomForcingField[2] == NULL)
	RandomForcingField[2] = new float[size];
      velocity3 = RandomForcingField[2];
    }
  else {
    velocity3 = new float[size];
    for (i = 0; i < size; i++)
      velocity3[i] = 0.0;
  }
  
  /* allows for the domain to start at position < 0 */
  
  for (dim = 0; dim < GridRank; dim++)
    Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
				(*(CellWidth[dim]))) - GridStartIndex[dim];
  
  /* compute linear box size in mesh points */

  nbox = nint((DomainRightEdge[0] - DomainLeftEdge[0])/
	      (*(CellWidth[0])));

   if (debug)
    printf("GCRFF: nbox = %d %d %d %d\n",
            nbox, Offset[0], Offset[1], Offset[2]);

  FORTRAN_NAME(turboinit)(&GridRank, &nbox, velocity1, velocity2, velocity3,
			  &GridDimension[0], &GridDimension[1], &GridDimension[2],
			  &Offset[0],&Offset[1],&Offset[2], &RandomForcingMachNumber);
                            
  /* deallocate temporary space for solver */
  
  if (GridRank < 3) delete [] velocity3;
  if (GridRank < 2) delete [] velocity2;
  
  return SUCCESS;

}


