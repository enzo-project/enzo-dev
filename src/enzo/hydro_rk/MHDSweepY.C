/***********************************************************************
/
/  COMPUTE 1D MHD FLUX IN Y DIRECTION
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

int MHDSweepY(float **Prim, float **Flux3D, int GridDimension[], 
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
    FluxLine[field] = new float[Yactivesize+1];
  }
  for (int field = 0; field < NEQ_MHD+NSpecies+NColor-idual; field++) {
    Prim1[field] = new float[GridDimension[1]];
  }

  int extra = (ReconstructionMethod == PPM);
  for (int field = 0; field < NEQ_MHD-idual; field++) {
    priml[field] = new float[Yactivesize+1+extra];
    primr[field] = new float[Yactivesize+1+extra];
  }

  for (int field = 0; field < NSpecies; field ++) {
    species[field] = new float[Yactivesize+1];
  }

  for (int field = 0; field < NColor; field ++) {
    colors[field] = new float[Yactivesize+1];
  }

  float min_coeff = 0.0;
  if (UseMinimumPressureSupport) {
    min_coeff = MinimumPressureSupportParameter*
      0.32*pow(CellWidth[1][0],2)/(Gamma*(Gamma-1.0));
  }

  float etot, vx, vy, vz, v2, p, Bx, By, Bz, B2, rho;
  for (k = 0; k < Zactivesize; k++) {
    for (i = 0; i < Xactivesize; i++) {

      // copy the relevant part of U and Prim into U1 and Prim1
      for (j = 0; j < GridDimension[1]; j++) {
	igrid = (i + GridStartIndex[0]) + j * GridDimension[0] +
	  (k + GridStartIndex[2]) * GridDimension[1] * GridDimension[0];
	
	rho = Prim[iden][igrid]; // density
	vx  = Prim[ivy ][igrid]; // vx = vy
	vy  = Prim[ivz ][igrid]; // vy = vz
	vz  = Prim[ivx ][igrid]; // vz = vx
	Bx  = Prim[iBy ][igrid];
	By  = Prim[iBz ][igrid];
	Bz  = Prim[iBx ][igrid];
	if (DualEnergyFormalism) {
	  Prim1[1][j] = Prim[ieint][igrid];
	} else {
	  etot = Prim[ietot][igrid];
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx*Bx + By*By + Bz*Bz;
	  Prim1[1][j] = etot - 0.5*v2 - 0.5*B2/rho;
	}

	if (EOSType > 0) {
	  float h, cs, dpdrho, dpde;
	  EOS(p, Prim[iden][igrid], Prim1[1][j], h, cs, dpdrho, dpde, EOSType, 0);
	  Prim1[1][j] = p;
	} 

	Prim1[1][j] = max(Prim1[1][j], min_coeff*rho);
	Prim1[0][j] = rho;
	Prim1[2][j] = vx;
	Prim1[3][j] = vy;
	Prim1[4][j] = vz;
        Prim1[5][j] = Bx;
        Prim1[6][j] = By;
        Prim1[7][j] = Bz;
        Prim1[8][j] = Prim[iPhi][igrid];
      }

      /* Copy species and color fields */

      for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++) {
	for (j = 0; j < GridDimension[1]; j++) {
	  igrid = (i + GridStartIndex[0]) + j * GridDimension[0] +
	    (k + GridStartIndex[2]) * GridDimension[1] * GridDimension[0];
	  Prim1[field-idual][j] = Prim[field][igrid];
	}
      }
	    
      // compute FluxLine from U1 and Prim1
      if (MHDLine(Prim1, priml, primr, species, colors, 
		  FluxLine, Yactivesize, dtdx, 'y', i, k, fallback) == FAIL) {
	printf("MHDLine failed.\n");
	return FAIL;
      }
      
      // copy FluxLine to the corresponding part of Flux3D
      for (j = 0; j < Yactivesize+1; j++) {
	iflux = i + (Xactivesize+1)*(j + k*(Yactivesize+1));
	Flux3D[iD   ][iflux] = FluxLine[iD   ][j];
	Flux3D[iS1  ][iflux] = FluxLine[iS3  ][j];
	Flux3D[iS2  ][iflux] = FluxLine[iS1  ][j];
	Flux3D[iS3  ][iflux] = FluxLine[iS2  ][j];
	Flux3D[iEtot][iflux] = FluxLine[iEtot][j];
	if (DualEnergyFormalism) {
	  Flux3D[iEint][iflux] = FluxLine[iEint][j];
	}
	Flux3D[iBx ][iflux] = FluxLine[iBz ][j];
	Flux3D[iBy ][iflux] = FluxLine[iBx ][j];
	Flux3D[iBz ][iflux] = FluxLine[iBy ][j];
	Flux3D[iPhi][iflux] = FluxLine[iPhi][j];
	for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++) {
	  Flux3D[field][iflux] = FluxLine[field][j];
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
