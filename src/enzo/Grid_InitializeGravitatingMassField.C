/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:  Daniel R. Reynolds 
/  date:       February, 2006
/
/  PURPOSE:
/
/  NOTE:
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
 
/* function prototypes */
 
 
int grid::InitializeGravitatingMassField(int RefinementFactor)
{
 
  /* Error check */
 
  if (RefinementFactor < 1 || RefinementFactor > RefineBy) {
    ENZO_VFAIL("RefinementFactor = %"ISYM" out of range.\n", RefinementFactor)
  }
 
  /* Check to see if the field was already initialized. */
 
  if (GravitatingMassFieldCellSize != FLOAT_UNDEFINED)
    return SUCCESS;
 
  int dim, GravityBufferSize = GRAVITY_BUFFER_SIZE, DimTemp, BufferSize;
 
  if( ! SelfGravity )
    GravityBufferSize = 0;
 
  /* Determine the size of the mass grid we'll need.
     1) For the top grid, this is just the active grid size (periodic)
     2) For the top grid, this is just twice the active grid size (isolated)
     3) For the subgrid we will use the boundary zones as well as the
        active region and then add some padding. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Make the GravitatingMassField the size of the active region
       plus the GravityBufferSize (in Parent cell units) on either size. */
 
    DimTemp = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    //      BufferSize = min(RefinementFactor*GravityBufferSize, DimTemp);
    BufferSize = RefinementFactor*GravityBufferSize;
    //      if (int(DimTemp/4)*4 != DimTemp && RefinementFactor == 2)

    BufferSize = ( (BufferSize <= NumberOfGhostZones ) ? NumberOfGhostZones + 1 : BufferSize ) ;
 
    GravitatingMassFieldDimension[dim] = DimTemp +
      2*max(BufferSize, NumberOfGhostZones);
    GravitatingMassFieldCellSize = CellWidth[dim][0];
    GravitatingMassFieldLeftEdge[dim] = GridLeftEdge[dim] -
      max(BufferSize, NumberOfGhostZones)*
      GravitatingMassFieldCellSize;
  }
 
  /* Set unused dims. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    GravitatingMassFieldDimension[dim] = 1;
    GravitatingMassFieldLeftEdge[dim] = DomainLeftEdge[dim];
  }


  return SUCCESS;
}
