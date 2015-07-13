/***********************************************************************
/
/  GRID CLASS (COMPUTE 3D FLUXES)
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
int HydroSweepX(float **Prim,  float **Flux3D, int GridDimension[], 
		int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback);
int HydroSweepY(float **Prim,  float **Flux3D, int GridDimension[], 
		int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback);
int HydroSweepZ(float **Prim,  float **Flux3D, int GridDimension[], 
		int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback);

int grid::Hydro3D(float **Prim, float **dU, float dt,
		  fluxes *SubgridFluxes[], int NumberOfSubgrids, 
		  float fluxcoef, int fallback)
  /* 
     Input:  U[NEQ_SRHYDRO][GridDimension^3]: the conserved variables vector 
                                              including ghost zones.
             Prim[NEQ_SRHYDRO+1][GridDimension^3]:
     Output: dU[NEQ_SRHYDRO][activesize^3]: the spatial differenced value without 
                                            ghost zones.
             SubgridFluxes:
  */
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int iflux, ifluxp1, fluxsize = 1, activesize = 1;
  float *Flux3D[NEQ_HYDRO+NSpecies+NColor];

  for (int dim = 0; dim < GridRank; dim++) {
    fluxsize *= (GridDimension[dim]-2*NumberOfGhostZones+1);
  }

  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
  }

  int Xactivesize = GridDimension[0]-2*NumberOfGhostZones;
  int Yactivesize = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  int Zactivesize = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    Flux3D[field] = new float[fluxsize];
  }
  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    for (int i = 0; i < fluxsize; i++) {
      Flux3D[field][i] = 0.0;
    }
  }

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

  // compute flux at cell faces in x direction
  float dtdx = dt/(a*CellWidth[0][0]);
  if (HydroSweepX(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback) == FAIL) {
    return FAIL;
  }

  // Update dU
  if (Coordinate == Cartesian) {
    for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
      int n = 0;
      for (int k = 0; k < Zactivesize; k++) {
	for (int j = 0; j < Yactivesize; j++) {
	  for (int i = 0; i < Xactivesize; i++, n++) {
	    iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
	    dU[field][n] = -(Flux3D[field][iflux+1]-Flux3D[field][iflux])*dtdx;
	  }
	}
      }
    }
  }

  FLOAT geofacr, geofacl;

  if (Coordinate == Spherical) {
    FLOAT xl, xr, xc;
    for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
      int n = 0, i1;
      for (int k = 0; k < Zactivesize; k++) {
        for (int j = 0; j < Yactivesize; j++) {
          for (int i = 0; i < Xactivesize; i++, n++) {
            iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
            i1 = i + GridStartIndex[0];
            xl = CellLeftEdge[0][i1];
            xc = xl + 0.5*CellWidth[0][i1];
            xr = xl + CellWidth[0][i1];
            geofacr = xr*xr/(xc*xc);
            geofacl = xl*xl/(xc*xc);
            dU[field][n] = -(geofacr*Flux3D[field][iflux+1]-geofacl*Flux3D[field][iflux])*dtdx;
          }
        }
      }
    }
  }


  if (FluxCorrection) {
    if (this->SaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				Flux3D, 0, fluxcoef, dt) == FAIL) {
      return FAIL;
    }
  }

  // Sweep Y

  if (GridRank > 1) {
    dtdx = dt/(a*CellWidth[1][0]);
    // compute flux in y direction
    if (HydroSweepY(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback) == FAIL)
      return FAIL;

    // Update dU
    if (Coordinate == Cartesian) {
      for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
	int n = 0;
	for (int k = 0; k < Zactivesize; k++) {
	  for (int j = 0; j < Yactivesize; j++) {
	    for (int i = 0; i < Xactivesize; i++, n++) {
	      iflux = i + (Xactivesize + 1) * (j+k*(Yactivesize + 1));
	      ifluxp1 = i + (Xactivesize + 1)*(j+1+k*(Yactivesize + 1));
	      dU[field][n] -= (Flux3D[field][ifluxp1]-Flux3D[field][iflux])*dtdx;
	    }
	  }
	}
      }

    }

    if (Coordinate == Spherical) {
      FLOAT yl, yr, yc, xc, sinyc;
      for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
        int n = 0, i1, j1;
        for (int k = 0; k < Zactivesize; k++) {
          for (int j = 0; j < Yactivesize; j++) {
            for (int i = 0; i < Xactivesize; i++, n++) {
              iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
              ifluxp1 = i + (Xactivesize + 1)*(j+1+k*(Yactivesize + 1));
              i1 = i + GridStartIndex[0];
              j1 = j + GridStartIndex[1];
              yl = CellLeftEdge[1][j1];
              yc = yl + 0.5*CellWidth[1][j1];
              yr = yl + CellWidth[1][j1];
              xc = CellLeftEdge[0][i1] + 0.5*CellWidth[0][i1];
              sinyc = sin(yc);
              geofacr = sin(yr)/(xc*sinyc);
              geofacl = sin(yl)/(xc*sinyc);
              dU[field][n] -= (geofacr*Flux3D[field][ifluxp1]-geofacl*Flux3D[field][iflux])*dtdx;
            }
          }
        }
      }
    }


    if (FluxCorrection) {
      if (this->SaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				  Flux3D, 1, fluxcoef, dt) == FAIL) {
	return FAIL;
      }
    }
  }
  
  // Sweep Z

  if (GridRank > 2) {
    dtdx = dt/(a*CellWidth[2][0]);
    // compute flux in z direction
    if (HydroSweepZ(Prim, Flux3D, GridDimension, GridStartIndex, CellWidth, dtdx, fallback) == FAIL)
      return FAIL;

    // Update dU
    if (Coordinate == Cartesian) {
      for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
	int n = 0;
	for (int k = 0; k < Zactivesize; k++) {
	  for (int j = 0; j < Yactivesize; j++) { 
	    for (int i = 0; i < Xactivesize; i++, n++) { 
	      iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
	      ifluxp1 = i + (Xactivesize + 1)*(j+(k + 1)*(Yactivesize + 1));
	      dU[field][n] -= (Flux3D[field][ifluxp1]-Flux3D[field][iflux])*dtdx;
	    }
	  }
	}
      }
    }

    if (Coordinate == Spherical) {
      FLOAT yc, xc, sinyc;
      for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
        int n = 0, i1, j1;
        for (int k = 0; k < Zactivesize; k++) {
          for (int j = 0; j < Yactivesize; j++) {
            for (int i = 0; i < Xactivesize; i++, n++) {
              iflux = i + (Xactivesize+1) * (j + k*(Yactivesize+1));
              ifluxp1 = i + (Xactivesize + 1)*(j+(k + 1)*(Yactivesize + 1));
              i1 = i + GridStartIndex[0];
              j1 = j + GridStartIndex[1];
              xc = CellLeftEdge[0][i1] + 0.5*CellWidth[0][i1];
              yc = CellLeftEdge[1][j1] + 0.5*CellWidth[1][j1];
              dU[field][n] -= (Flux3D[field][ifluxp1] - Flux3D[field][iflux])*dtdx/(xc*sin(yc));
            }
          }
        }
      }
    }


    if (FluxCorrection) {
      if (this->SaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids,
				  Flux3D, 2, fluxcoef, dt) == FAIL) {
	return FAIL;
      }
    }
  }


  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    delete [] Flux3D[field];
  }

  return SUCCESS;

}
