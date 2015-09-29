/***********************************************************************
/
/  COMPUTE FLUX IN Z DIRECTION
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
#include "EOS.h"
#include "phys_constants.h"

int HydroLine(float **Prim, float **priml, float **primr,
	      float **species, float **colors, float **FluxLine, int ActiveSize,
	      float dtdx, char direc, int ij, int ik, int fallback);

int HydroSweepZ(float **Prim, float **Flux3D, int GridDimension[], 
		int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback)
  /*
    Input: U[NEQ_HYDRO][GridDimension^3].
           Prim[NEQ_HYDRO+1][GridDimension^3].
    Output: Flux3D[NEQ_HYDRO][(activesize+1)^3]
  */
{

  int i, j, k, m, iflux, igrid;
  int idual = (DualEnergyFormalism) ? 1 : 0;
  float *FluxLine[NEQ_HYDRO+NSpecies+NColor];
  float *Prim1[NEQ_HYDRO+NSpecies+NColor-idual];
  float *priml[NEQ_HYDRO-idual], *primr[NEQ_HYDRO-idual], *species[NSpecies], *colors[NColor];
  
  int Xactivesize = GridDimension[0]-2*NumberOfGhostZones;
  int Yactivesize = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  int Zactivesize = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    FluxLine[field] = new float[Zactivesize+1];
  }

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor-idual; field++) {
    Prim1[field] = new float[GridDimension[2]];
  }

  int extra = (ReconstructionMethod == PPM);
  for (int field = 0; field < NEQ_HYDRO-idual; field++) {
    priml[field] = new float[Zactivesize+1+extra];
    primr[field] = new float[Zactivesize+1+extra];
  }

  for (int field = 0; field < NSpecies; field ++) {
    species[field] = new float[Zactivesize+1];
  }

  for (int field = 0; field < NColor; field ++) {
    colors[field] = new float[Zactivesize+1];
  }


 
  float min_coeff = 0.0;
  if (UseMinimumPressureSupport) {
    min_coeff = GravitationalConstant/(4.0*pi) / (pi * (Gamma*(Gamma-1.0))) *
      MinimumPressureSupportParameter *
      CellWidth[0][0] * CellWidth[0][0];
  }

  float etot, vx, vy, vz, v2, p;
  for (j = 0; j < Yactivesize; j++) {
    for (i = 0; i < Xactivesize; i++) {

      // copy the relevant part of U and Prim into U1 and Prim1      
      for (k = 0; k < GridDimension[2]; k++) {
	igrid = (i + GridStartIndex[0]) + (j+GridStartIndex[1]) * GridDimension[0] +
	  k * GridDimension[1] * GridDimension[0];
	Prim1[0][k] = Prim[iden ][igrid]; // density
	vx = Prim[ivz  ][igrid]; // vx = vz
	vy = Prim[ivx  ][igrid]; // vy = vx
	vz = Prim[ivy  ][igrid]; // vz = vy

	if (DualEnergyFormalism) {
	  Prim1[1][k] = Prim[ieint][igrid];
	}
	else {
	  etot = Prim[ietot][igrid];
	  v2 = vx*vx + vy*vy + vz*vz;
	  Prim1[1][k] = etot - 0.5*v2;
	}
	if (EOSType > 0) {
	  float h, cs, dpdrho, dpde;
	  EOS(p, Prim[iden][igrid], Prim1[1][k], h, cs, dpdrho, dpde, EOSType, 0);
	  Prim1[1][k] = p;
	  // then compare pressures, not energies, if using floor
	  Prim1[1][k] = max(Prim1[1][k], min_coeff*Prim1[0][k]*Prim1[0][k]*(Gamma-1.0));
	}
	else
	  // compare energies if using floor
	  Prim1[1][k] = max(Prim1[1][k], min_coeff*Prim1[0][k]);

	Prim1[2][k] = vx;
	Prim1[3][k] = vy;
	Prim1[4][k] = vz;
      }

      for (int field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++) {
	for (k = 0; k < GridDimension[2]; k++) {
	  igrid = (i + GridStartIndex[0]) + (j+GridStartIndex[1]) * GridDimension[0] +
	    k * GridDimension[1] * GridDimension[0];
	  Prim1[field-idual][k] = Prim[field][igrid];
	}
      }

      // compute FluxLine from U1 and Prim1
      if (HydroLine(Prim1, priml, primr, species, colors, 
		    FluxLine, Zactivesize, dtdx, 'z', i, j, fallback) == FAIL) {
	printf("HydroLine failed in SweepZ\n");
	return FAIL;
      }

      // copy FluxLine to the corresponding part of Flux3D
      for (k = 0; k < Zactivesize+1; k++) {
	iflux = i + (Xactivesize+1)*(j + k*(Yactivesize+1));
	Flux3D[iD  ][iflux] = FluxLine[iD  ][k];
	Flux3D[iS1 ][iflux] = FluxLine[iS2 ][k];
	Flux3D[iS2 ][iflux] = FluxLine[iS3 ][k];
	Flux3D[iS3 ][iflux] = FluxLine[iS1 ][k];
	Flux3D[iEtot][iflux] = FluxLine[iEtot][k];
	if (DualEnergyFormalism) {
	  Flux3D[iEint][iflux] = FluxLine[iEint][k];
	}
	for (int field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++) {
	  Flux3D[field][iflux] = FluxLine[field][k];
	}
      }
    }
  }

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    delete [] FluxLine[field];
  }

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor-idual; field++) {
    delete [] Prim1[field];
  }

  for (int field = 0; field < NEQ_HYDRO-idual; field++) {
    delete [] priml[field];
    delete [] primr[field];
  }

  for (int field = 0; field < NSpecies; field++) {
    delete [] species[field];
  }

  for (int field = 0; field < NColor; field++) {
    delete [] colors[field];
  }


  return SUCCESS;
}
