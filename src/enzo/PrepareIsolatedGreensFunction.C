/***********************************************************************
/
/  PREPARES THE GREENS FUNCTION IN REAL SPACE
/
/  written by: Greg Bryan
/  date:       March, 2001
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

//  Compute derived quantites
//

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
#include "TopGridData.h"


int PrepareIsolatedGreensFunction(region *GreensFunction, int proc, 
				  int DomainDim[], TopGridData *MetaData)
{

  /* Declarations. */

  int i, j, k, dim, GridRank = MetaData->TopGridRank;
  float x, DomainCellSize;

  /* Create DomainDim and double it. */

  for (dim = 0; dim < GridRank; dim++)
    DomainDim[dim] = MetaData->TopGridDims[dim]*2;
  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    DomainDim[dim] = 1;

  /* Add two to the x-direction as required by the real-to-complex FFT. */

  DomainDim[0] += 2;

  /* Compute the domain cell size. */

  DomainCellSize = 1.0/float(DomainDim[0]);

  /* Create a slab of size DomainDim[0]/nproc * DomainDim[1] * DomainDim[2]. */

  x = float(proc) / float(NumberOfProcessors);
  GreensFunction->StartIndex[0] = nint(x/(2.0*DomainCellSize))*2;
  x += 1.0/float(NumberOfProcessors);
  GreensFunction->RegionDim[0] = nint(x/(2.0*DomainCellSize))*2 - 
                                 GreensFunction->StartIndex[0];
  for (dim = 1; dim < MAX_DIMENSION; dim++) {
    GreensFunction->StartIndex[dim] = 0;
    GreensFunction->RegionDim[dim] = DomainDim[dim];
  }

  GreensFunction->Processor = proc;
  GreensFunction->Data = NULL;

  /* Return if this is not our processor. */

  //  fprintf(stderr, "%"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", proc, MyProcessorNumber, GreensFunction->RegionDim[0], GreensFunction->RegionDim[1], GreensFunction->RegionDim[2]);
  if (proc != MyProcessorNumber)
    return SUCCESS;

  /* Compute size and allocate field with size of GravitatingMassField. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GreensFunction->RegionDim[dim];

  GreensFunction->Data = new float[size];

  /* Compute the real cell widths. */

  float RealCellWidth[MAX_DIMENSION], RealCellVolume = 1.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    RealCellWidth[dim] = (DomainRightEdge[dim]-DomainLeftEdge[dim])/
      float(MetaData->TopGridDims[dim]);
    if (dim < GridRank)
      RealCellVolume *= RealCellWidth[dim];
  }

  /* Set the constant to be used. */

  float GravConst, pi = M_PI;
  if (GridRank == 3) 
    GravConst = -GravitationalConstant*RealCellVolume/(4.0*pi);
  if (GridRank == 2)
    GravConst = GravitationalConstant*RealCellVolume/(2.0*pi);
  if (GridRank == 1)
    GravConst = GravitationalConstant*RealCellVolume/2.0;

  /* Set Greens' function. */

  int n = 0, i1, j1, k1;
  float xpos, ypos, zpos, r;
  for (k = 0; k < GreensFunction->RegionDim[2]; k++) {
    k1 = k+GreensFunction->StartIndex[2];
    if (k1 >= MetaData->TopGridDims[2])
      k1 = k1 - 2*MetaData->TopGridDims[2];
    zpos = (GridRank > 2) ? ((float(k1))*RealCellWidth[2]) : 0;
    for (j = 0; j < GreensFunction->RegionDim[1]; j++) {
      j1 = j+GreensFunction->StartIndex[1];
      if (j1 >= MetaData->TopGridDims[1])
	j1 = j1 - 2*MetaData->TopGridDims[1];
      ypos = (GridRank > 1) ? ((float(j1))*RealCellWidth[1]) : 0;
      for (i = 0; i < GreensFunction->RegionDim[0]; i++, n++) {
	i1 = i+GreensFunction->StartIndex[0];
	if (i1 >= MetaData->TopGridDims[0])
	  i1 = i1 - 2*MetaData->TopGridDims[0];
	xpos = (float(i1))*RealCellWidth[0];
	r = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
	r = max(r, 0.38*RealCellWidth[0]);
	//	r *= GravitatingMassFieldCellSize;
	if (GridRank == 3)
	  GreensFunction->Data[n] = GravConst/r;
	if (GridRank == 2)
	  GreensFunction->Data[n] = GravConst*log(r);
	if (GridRank == 1)
	  GreensFunction->Data[n] = GravConst*r;
      }
    }
  }



  return SUCCESS;
}
