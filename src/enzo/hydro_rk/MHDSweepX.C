/***********************************************************************
/
/  COMPUTE 1D MHD FLUX IN X DIRECTION
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

int MHDSweepX(float **Prim, float **Flux3D, int GridDimension[], 
	      int GridStartIndex[], FLOAT **CellWidth, float dtdx, int fallback)
  /*
    Input: U[NEQ_MHD][GridDimension^3].
           Prim[NEQ_MHD+1][GridDimension^3].
    Output: Flux3D[NEQ_MHD][(activesize+1)^3]
  */
{

  int i, j, k, m, iflux, igrid;
  int idual = (DualEnergyFormalism) ? 1 : 0;
//  int icons = (ConservativeReconstruction) ? 1 : 0;  // not implemented properly yet, TA
  float *FluxLine[NEQ_MHD+NSpecies+NColor];
  float *Prim1[NEQ_MHD+NSpecies+NColor-idual]; 
  float *priml[NEQ_MHD-idual], *primr[NEQ_MHD-idual], *species[NSpecies], *colors[NColor];
  
  int Xactivesize = GridDimension[0]-2*NumberOfGhostZones;
  int Yactivesize = GridDimension[1] > 1 ? GridDimension[1]-2*NumberOfGhostZones : 1;
  int Zactivesize = GridDimension[2] > 1 ? GridDimension[2]-2*NumberOfGhostZones : 1;


  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    FluxLine[field] = new float[Xactivesize+1];
  }

  //+icons may be needed for ConservativeReconstruction implementation
  for (int field = 0; field < NEQ_MHD+NSpecies+NColor+idual; field++) {
    Prim1[field] = new float[GridDimension[0]];
  }

  int extra = (ReconstructionMethod == PPM);
  for (int field = 0; field < NEQ_MHD-idual; field++) {
    priml[field] = new float[Xactivesize+1+extra];
    primr[field] = new float[Xactivesize+1+extra];
  }

  for (int field = 0; field < NSpecies; field ++) {
    species[field] = new float[Xactivesize+1];
  }

  for (int field = 0; field < NColor; field ++) {
    colors[field] = new float[Xactivesize+1];
  }

  float etot, vx, vy, vz, v2, p, Bx, By, Bz, B2, rho;
  
  float min_coeff = 0.0;
  if (UseMinimumPressureSupport) {
    min_coeff = MinimumPressureSupportParameter*
      0.32*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));
  }


  for (k = 0; k < Zactivesize; k++) {
    for (j = 0; j < Yactivesize; j++) {

      for (i = 0; i < GridDimension[0]; i++) {
	// We only use density, internal energy, vi to do the reconstruction
	igrid = i + (j + GridStartIndex[1]
		     + (k + GridStartIndex[2])*GridDimension[1]) * GridDimension[0];
	rho = Prim[iden][igrid];
	vx  = Prim[ivx ][igrid];
	vy  = Prim[ivy ][igrid];
	vz  = Prim[ivz ][igrid];
	Bx  = Prim[iBx ][igrid];
	By  = Prim[iBy ][igrid];
	Bz  = Prim[iBz ][igrid];

	if (DualEnergyFormalism) {
	  Prim1[1][i] = Prim[ieint][igrid];
	}
	else {
	  etot = Prim[ietot][igrid];
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx*Bx + By*By + Bz*Bz;
	  Prim1[1][i] = etot - 0.5*v2 - 0.5*B2/rho;
	}

	if (EOSType > 0) {
	  float h, cs, dpdrho, dpde;
	  EOS(p, Prim[iden][igrid], Prim1[1][i], h, cs, dpdrho, dpde, EOSType, 0);
	  Prim1[1][i] = p;
	}

	Prim1[0][i] = rho;
	Prim1[1][i] = max(Prim1[1][i], min_coeff*rho);
	Prim1[2][i] = vx;
	Prim1[3][i] = vy;
	Prim1[4][i] = vz;
	Prim1[5][i] = Bx;
	Prim1[6][i] = By;
	Prim1[7][i] = Bz;
	Prim1[8][i] = Prim[iPhi][igrid];
//	if (ConservativeReconstruction)  // add dx/dt for every cell
//	  Prim1[9][i] = dtdx*CellWidth[0][0]/CellWidth[0][i];
      }

      /* Copy species and color fields */
      
      for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++) {
	for (i = 0; i < GridDimension[0]; i++) {
	  igrid = i + (j + GridStartIndex[1]
		       + (k + GridStartIndex[2])*GridDimension[1]) * GridDimension[0];
	  Prim1[field-idual][i] = Prim[field][igrid];
	}
      }

      // compute FluxLine from U1 and Prim1
      if (MHDLine(Prim1, priml, primr, species, colors, 
		  FluxLine, Xactivesize, dtdx, 'x', j, k, fallback) == FAIL) {
	printf("MHDSweepX: MHDLine failed.\n");
	return FAIL;
      }

      // copy FluxLine to the corresponding part of Flux3D
      for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
	for (i = 0; i < Xactivesize+1; i++) {
	  iflux = i + (Xactivesize+1)*(j + k*(Yactivesize+1));
	  Flux3D[field][iflux] = FluxLine[field][i];
	}
      }
    }
  }

	     

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] FluxLine[field];
  }
  //+icons may be needed for ConservativeReconstruction implementation
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

