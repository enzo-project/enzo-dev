/***********************************************************************
/
/  GRID CLASS (SAVE SUBGRID BOUNDARY FLUXES)
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::SaveSubgridFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids,
			    float *Flux3D[], int flux, float fluxcoef, float dt)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int Start[MAX_DIMENSION], End[MAX_DIMENSION], Dim[MAX_DIMENSION], Activesize[MAX_DIMENSION],
    Offset, iflux, igridflux;


  Activesize[0] = GridDimension[0]-2*NumberOfGhostZones;
  Activesize[1] = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  Activesize[2] = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;

  FLOAT a = 1, dadt;

  /* If using comoving coordinates, multiply dx by a(n+1/2).
     In one fell swoop, this recasts the equations solved by solver
     in comoving form (except for the expansion terms which are taken
     care of elsewhere). */
  
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
      return FAIL;
    }
  }
  
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

    //printf("Start=%"ISYM", End=%"ISYM"\n", Start[flux], End[flux]);

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      Dim[dim] = End[dim] - Start[dim] + 1;
    }

    Offset = SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[flux][flux]
      - SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][flux] + 1;
    //Offset = min(Offset, GridDimension[flux] - 1);

    for (int dim = 0; dim < flux; dim++) {
      Offset *= (Activesize[dim]+1);
    }
    
    FLOAT dtdx = dtFixed/(a*CellWidth[flux][0]);
    for (int field = 0; field < NEQ_HYDRO; field++) {
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
	    //printf("iflux=%"ISYM", igridflux=%"ISYM", offset=%"ISYM"\n",
	    //   iflux, igridflux, Offset);
	  }
	}
      }
    }
  }
    
  /*for (int field = 0; field < NEQ_SRHYDRO; field++) {
    for (int n = 0; n < Activesize[0]+1; n++) {
      printf("%"GSYM" ", Flux3D[field][n]);
    }
    printf("\n");
    }*/
  /*printf("EstimateFluxes: left=%"GSYM", right=%"GSYM"\n",
	 SubgridFluxes[0]->LeftFluxes[1][0][0],
	 SubgridFluxes[0]->RightFluxes[1][0][0]);*/

  return SUCCESS;

}
