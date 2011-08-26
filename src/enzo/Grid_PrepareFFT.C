/***********************************************************************
/
/  PREPARE FOR A PARALLEL FFT
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
 
 
int grid::PrepareFFT(region *InitialRegion, int Field, int DomainDim[])
{
 
  int dim, size;
  float DomainCellSize[MAX_DIMENSION];
 
  /* Set Domain parameters. */
 
  for (dim = 0; dim < GridRank; dim++) {
    DomainCellSize[dim] = GravitatingMassFieldCellSize;
    DomainDim[dim] = nint((DomainRightEdge[dim] - DomainLeftEdge[dim])/
                          DomainCellSize[dim]);
  }
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    DomainCellSize[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    DomainDim[dim] = 1;
  }
 
  /* -------------------------------------------------- */
  /* Generate a description of the initial data layout. */
 
  int GravDim[] = {1,1,1}, Zero[] = {0,0,0}, GravStart[] = {0,0,0};
 
  /* Set start index and dimension of region,
     and add 2 if we are on the right-hand side of the x-dimension,
     which is required for the real-to-complex FFT. */
 
  for (dim = 0, size = 1; dim < GridRank; dim++) {
    InitialRegion->StartIndex[dim] = nint((GridLeftEdge[dim] -
	       DomainLeftEdge[dim])/DomainCellSize[dim]);
    GravDim[dim] = GravitatingMassFieldDimension[dim];
    GravStart[dim] = nint((GridLeftEdge[dim] -
	       GravitatingMassFieldLeftEdge[dim])/DomainCellSize[dim]);
    InitialRegion->RegionDim[dim] = GridEndIndex[dim] -
	                               GridStartIndex[dim] + 1;
    if (InitialRegion->StartIndex[dim] +
	InitialRegion->RegionDim[dim] == DomainDim[dim] &&
	dim == 0) {
      InitialRegion->RegionDim[dim] += 2;
      //      InitialRegion->RegionDim[dim] += 0;
    }
    size *= InitialRegion->RegionDim[dim];
  }
  DomainDim[0] += 2;
 
  /* set start index and region for unused dims. */
 
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    InitialRegion->StartIndex[dim] = 0;
    InitialRegion->RegionDim[dim]  = 1;
  }
 
  /* Set processor number where data lives. */
 
  InitialRegion->Processor = ProcessorNumber;
  InitialRegion->Data      = NULL;
 
  /* If the data is on this processor then copy it to a new region. */
 
  if (MyProcessorNumber == InitialRegion->Processor) {
 
    /* Set FieldPointer to the appropriate field. */
 
    float *FieldPointer = NULL;
    if (Field == GRAVITATING_MASS_FIELD)
      FieldPointer = GravitatingMassField;
    if (Field == POTENTIAL_FIELD)
      FieldPointer = PotentialField;
    if (FieldPointer == NULL) {
      ENZO_VFAIL("Field type %"ISYM" not recognized.\n", Field)
    }
 
    InitialRegion->Data = new float[size];
 
    FORTRAN_NAME(copy3d)(FieldPointer, InitialRegion->Data,
			 GravDim, GravDim+1, GravDim+2,
			 InitialRegion->RegionDim,
			 InitialRegion->RegionDim+1,
			 InitialRegion->RegionDim+2,
			 Zero, Zero+1, Zero+2,
			 GravStart, GravStart+1, GravStart+2);
 
    /* Delete old field. */
 
    if (Field == GRAVITATING_MASS_FIELD) {
      //      delete FieldPointer;
      //      GravitatingMassField = NULL;
    }
    if (Field == POTENTIAL_FIELD) {
      delete FieldPointer;
      PotentialField = NULL;
    }
 
  } // end: if (MyProcessorNumber == ...)

 
  return SUCCESS;
}
