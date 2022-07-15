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

  if (GravitySolverType == GRAVITY_SOLVER_FAST) {
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
  } else if (GravitySolverType == GRAVITY_SOLVER_APM) {

    /* Determine the size of the mass grid we'll need.
       - For the top grid, this is just the active grid size plus a buffer
       - For the subgrid we will use the boundary zones as well as the
       active region and then add some padding. */

    for (dim = 0; dim < GridRank; dim++) {
      switch (GravityBoundaryType) {

      /* 1) TopGrid Periodic or Isolated */
      case TopGridPeriodic:
      case TopGridIsolated:
        DimTemp = GridEndIndex[dim] - GridStartIndex[dim] + 1;
        BufferSize = RefinementFactor*GravityBufferSize;
        BufferSize = ( (BufferSize <= NumberOfGhostZones ) ? NumberOfGhostZones + 1 : BufferSize ) ;

        GravitatingMassFieldDimension[dim] = DimTemp +
          2*max(BufferSize, NumberOfGhostZones);
        GravitatingMassFieldCellSize = CellWidth[dim][0];
        GravitatingMassFieldLeftEdge[dim] = GridLeftEdge[dim] -
          max(BufferSize, NumberOfGhostZones)*
          GravitatingMassFieldCellSize;
        break;

      /* 3) Subgrid */
      case SubGridIsolated:{

        /* Compute the extra padding required to include all the mass
           within one convolution kernal radius of the cells on the edge.
           This is some fraction of parent grid's particle smoothing size
           minues whatever buffer is already there. */

        int SubGridExtra = max(nint(float(RefinementFactor)*S2ParticleSize*0.65 -
                                    float(GravityResolution)*
                                    float(GridStartIndex[dim]) ), 0);
        GravitatingMassFieldDimension[dim] =
          nint(float(GridDimension[dim])*GravityResolution)  +
          //              nint(float(RefinementFactor)*S2ParticleSize)       +
          //              FFT_SAFETY_FACTOR                                  +
          2*SubGridExtra;
        GravitatingMassFieldCellSize = CellWidth[dim][0]/GravityResolution;
        GravitatingMassFieldLeftEdge[dim] = CellLeftEdge[dim][0] -
          float(SubGridExtra)*GravitatingMassFieldCellSize;
        break;
      }

      /* 4) undefined or unknown is an error */
      case GravityUndefined:
      default:
        fprintf(stderr, "GravityBoundaryType undefined.\n");
        return FAIL;

      } // end switch
    } // end loop over dims
  } // end: if (GravitySolverType == GRAVITY_SOLVER_FAST)
 
  /* Set unused dims. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    GravitatingMassFieldDimension[dim] = 1;
    GravitatingMassFieldLeftEdge[dim] = DomainLeftEdge[dim];
  }


  return SUCCESS;
}
