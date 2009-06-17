/***********************************************************************
/
/  GRID CLASS (SAVE SUBGRID BOUNDARY FLUXES)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::SaveMHDSubgridFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids,
			       float *Flux3D[], int flux, float fluxcoef, float dt)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int Start[MAX_DIMENSION], End[MAX_DIMENSION], Dim[MAX_DIMENSION], Activesize[MAX_DIMENSION],
    Offset, iflux, igridflux;


  Activesize[0] = GridDimension[0]-2*DEFAULT_GHOST_ZONES;
  Activesize[1] = GridDimension[1] > 1 ? GridDimension[1]-2*DEFAULT_GHOST_ZONES : 1;
  Activesize[2] = GridDimension[2] > 1 ? GridDimension[2]-2*DEFAULT_GHOST_ZONES : 1;

  
  /* Start: subgrid left flux start index relative to this grid
     End: subgrid left flux start index relative to this grid */

  for (int subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      Start[dim] = 0;
      End[dim] = 0;
    }

    for (int dim = 0; dim < GridRank; dim++) {
      Start[dim] = SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][dim]
	-nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/CellWidth[dim][0]);
      End[dim] = SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[flux][dim]
	-nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/CellWidth[dim][0]);

      //Start[flux] = max(Start[flux]-1, 0);
      End[flux] = Start[flux];

    }

    //printf("Start=%d, End=%d\n", Start[flux], End[flux]);

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      Dim[dim] = End[dim] - Start[dim] + 1;
    }

    Offset = SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[flux][flux]
      - SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][flux] + 1;
    //Offset = min(Offset, GridDimension[flux] - 1);

    for (int dim = 0; dim < flux; dim++) {
      Offset *= (Activesize[dim]+1);
    }
    
    FLOAT dtdx = dtFixed/CellWidth[flux][0];
    for (int field = 0; field < NEQ_MHD; field++) {
      for (int k = Start[2]; k <= End[2]; k++) {
	for (int j = Start[1]; j <= End[1]; j++) {
	  for (int i = Start[0]; i <= End[0]; i++) {
	    iflux = (k-Start[2])*Dim[1]*Dim[0] + (j-Start[1])*Dim[0] +
	      (i-Start[0]);
	    igridflux = i + j * (Activesize[0] + 1) +
	      k * (Activesize[0] + 1) * (Activesize[1] + 1);
	    SubgridFluxes[subgrid]->LeftFluxes[field][flux][iflux] += 
	      fluxcoef*Flux3D[field][igridflux]*dtdx;
	    SubgridFluxes[subgrid]->RightFluxes[field][flux][iflux] += 
	      fluxcoef*Flux3D[field][igridflux+Offset]*dtdx;
	    //printf("iflux=%d, igridflux=%d, offset=%d\n",
	    //   iflux, igridflux, Offset);
	  }
	}
      }
    }
  }

    
  /*for (int field = 0; field < NEQ_SRHYDRO; field++) {
    for (int n = 0; n < Activesize[0]+1; n++) {
      printf("%g ", Flux3D[field][n]);
    }
    printf("\n");
    }*/
  /*printf("EstimateFluxes: left=%g, right=%g\n",
	 SubgridFluxes[0]->LeftFluxes[1][0][0],
	 SubgridFluxes[0]->RightFluxes[1][0][0]);*/


  return SUCCESS;

}
