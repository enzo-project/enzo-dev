/***********************************************************************
/
/  COMPUTE 1D MHD FLUX IN Z DIRECTION
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

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "EOS.h"

int MHDLine(float **Prim, float **priml, float **primr,
	    float **species, float **colors, float **FluxLine, int ActiveSize,
	    float dtdx, char direc, int jj, int kk, int fallback);

int MHDSweepZ(float **Prim, float **Flux3D, int GridDimension[], 
	      int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback)
  /*
    Input: U[NEQ_MHD][GridDimension^3].
           Prim[NEQ_MHD+1][GridDimension^3].
    Output: Flux3D[NEQ_MHD][(activesize+1)^3]
  */
{

  int i, j, k, m, iflux, igrid;
  int idual = (DualEnergyFormalism) ? 1 : 0;
  float *FluxLine[NEQ_MHD+NSpecies+NColor];
  float *Prim1[NEQ_MHD+NSpecies+NColor-idual];
  float *priml[NEQ_MHD-idual], *primr[NEQ_MHD-idual], *species[NSpecies], *colors[NColor];
  
  int Xactivesize = GridDimension[0]-2*NumberOfGhostZones;
  int Yactivesize = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  int Zactivesize = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    FluxLine[field] = new float[Zactivesize+1];
  }
  for (int field = 0; field < NEQ_MHD+NSpecies+NColor-idual; field++) {
    Prim1[field] = new float[GridDimension[2]];
  }

  int extra = (ReconstructionMethod == PPM);
  for (int field = 0; field < NEQ_MHD-idual; field++) {
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
    min_coeff = MinimumPressureSupportParameter*
      0.32*pow(CellWidth[2][0],2)/(Gamma*(Gamma-1.0));
  }

  float etot, vx, vy, vz, v2, p, rho, Bx, By, Bz, B2;
  for (j = 0; j < Yactivesize; j++) {
    for (i = 0; i < Xactivesize; i++) {

      // copy the relevant part of U and Prim into U1 and Prim1      
      for (k = 0; k < GridDimension[2]; k++) {
	igrid = (i + GridStartIndex[0]) + (j+GridStartIndex[1]) * GridDimension[0] +
	  k * GridDimension[1] * GridDimension[0];
	rho = Prim[iden][igrid]; // density
	vx  = Prim[ivz ][igrid]; // vx = vz
	vy  = Prim[ivx ][igrid]; // vy = vx
	vz  = Prim[ivy ][igrid]; // vz = vy
	Bx  = Prim[iBz ][igrid];
	By  = Prim[iBx ][igrid];
	Bz  = Prim[iBy ][igrid];

	if (DualEnergyFormalism) {
	  Prim1[1][k] = Prim[ieint][igrid];
	} else {
	  etot = Prim[ietot][igrid];
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx*Bx + By*By + Bz*Bz;
	  Prim1[1][k] = etot - 0.5*v2 - 0.5*B2/rho;
	}

	if (EOSType > 0) {
	  float h, cs, dpdrho, dpde;
	  EOS(p, Prim[iden][igrid], Prim1[1][k], h, cs, dpdrho, dpde, EOSType, 0);
	  Prim1[1][k] = p;
	} 

	Prim1[1][k] = max(Prim1[1][k], min_coeff*rho);
	Prim1[0][k] = rho;
	Prim1[2][k] = vx;
	Prim1[3][k] = vy;
	Prim1[4][k] = vz;
	Prim1[5][k] = Bx;
        Prim1[6][k] = By;
        Prim1[7][k] = Bz;
        Prim1[8][k] = Prim[iPhi][igrid];
      }

      /* Copy species and color fields */

      for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++) {
	for (k = 0; k < GridDimension[2]; k++) {
	  igrid = (i + GridStartIndex[0]) + (j+GridStartIndex[1]) * GridDimension[0] +
	    k * GridDimension[1] * GridDimension[0];
	  Prim1[field-idual][k] = Prim[field][igrid];
	}
      }

      // compute FluxLine from U1 and Prim1
      if (MHDLine(Prim1, priml, primr, species, colors, 
		  FluxLine, Zactivesize, dtdx, 'z', i, j, fallback) == FAIL) {
	printf("MHDLine failed in SweepZ\n");
	return FAIL;
      }

      // copy FluxLine to the corresponding part of Flux3D
      for (k = 0; k < Zactivesize+1; k++) {
	iflux = i + (Xactivesize+1)*(j + k*(Yactivesize+1));
	Flux3D[iD   ][iflux] = FluxLine[iD  ][k];
	Flux3D[iS1  ][iflux] = FluxLine[iS2 ][k];
	Flux3D[iS2  ][iflux] = FluxLine[iS3 ][k];
	Flux3D[iS3  ][iflux] = FluxLine[iS1 ][k];
	Flux3D[iEtot][iflux] = FluxLine[iEtot][k];
	if (DualEnergyFormalism) {
	  Flux3D[iEint][iflux] = FluxLine[iEint][k];
	}
	Flux3D[iBx ][iflux] = FluxLine[iBy][k];
	Flux3D[iBy ][iflux] = FluxLine[iBz][k];
	Flux3D[iBz ][iflux] = FluxLine[iBx][k];
	Flux3D[iPhi][iflux] = FluxLine[iPhi][k];
	for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++) {
	  Flux3D[field][iflux] = FluxLine[field][k];
	}
      }
    }
  }

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] FluxLine[field];
  }
  for (int field = 0; field < NEQ_MHD+NSpecies+NColor-idual; field++) {
    delete [] Prim1[field];
  }
  for (int field = 0; field < NEQ_MHD-idual; field++) {
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
