/***********************************************************************
/
/  WRAP-UP A PARALLEL FFT (COPY POTENTIAL BACK TO GRID)
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
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
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
 
int grid::FinishFFT(region *InitialRegion, int Field, int DomainDim[])
{
 
  int dim, size;
  float DomainCellSize[MAX_DIMENSION];
 
  /* Set Domain parameters. */
 
  for (dim = 0; dim < GridRank; dim++)
    DomainCellSize[dim] = GravitatingMassFieldCellSize;
  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    DomainCellSize[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
 
  /* -------------------------------------------------- */
  /* Generate a description of the initial data layout. */
 
  int GravDim[] = {1,1,1}, Zero[] = {0,0,0}, GravStart[] = {0,0,0};
 
  /* Set start index and dimension of region,
     and add 2 if we are on the right-hand side of the x-dimension,
     which is required for the real-to-complex FFT. */
 
  for (dim = 0, size = 1; dim < GridRank; dim++) {
    GravDim[dim] = GravitatingMassFieldDimension[dim];
    GravStart[dim] = nint((GridLeftEdge[dim] -
	       GravitatingMassFieldLeftEdge[dim])/DomainCellSize[dim]);
    size *= GravDim[dim];
  }
 
  /* If the data is on this processor then copy it to a new region. */
 
  if (MyProcessorNumber == InitialRegion->Processor) {
 
    /* Set FieldPointer to the appropriate field. */
 
    float *FieldPointer;
    if (Field == POTENTIAL_FIELD) {
      if (PotentialField == NULL)
	PotentialField = new float[size]();
      FieldPointer = PotentialField;
    } else {
      ENZO_VFAIL("Field %"ISYM" not recognized.\n", Field)
    }
 
    /* Copy region data into grid. */
 
    FORTRAN_NAME(copy3d)(InitialRegion->Data, FieldPointer,
			 InitialRegion->RegionDim,
			 InitialRegion->RegionDim+1,
			 InitialRegion->RegionDim+2,
			 GravDim, GravDim+1, GravDim+2,
			 GravStart, GravStart+1, GravStart+2,
			 Zero, Zero+1, Zero+2);
 
    /* Delete old field. */
 
    delete [] InitialRegion->Data;
    InitialRegion->Data = NULL;
 
  } // end: if (MyProcessorNumber == ...)

 
  return SUCCESS;
}
