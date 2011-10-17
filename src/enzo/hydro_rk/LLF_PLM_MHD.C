/***********************************************************************
/
/  LLF-PLM MHD SOLVER
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
#include "../hydro_rk/ReconstructionRoutines.h"

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq);
int plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize);
int plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize);
#ifdef ECUDA
void cuda_plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq, float Theta_Limiter);
int cuda_plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize, float Theta_Limiter, int NSpecies);
int cuda_plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize, float Theta_Limiter, int NColor);
#endif
int llf_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize);

int LLF_PLM_MHD(float **prim, float **priml, float **primr,
		float **species, float **colors,  float **FluxLine, int ActiveSize,
		char direc, int jj, int kk)
{

  // compute priml and primr
#ifdef ECUDA
  if (UseCUDA) {
    cuda_plm(prim, priml, primr, ActiveSize, 9, Theta_Limiter);
    } else   
#endif
  if (plm(prim, priml, primr, ActiveSize, 9) == FAIL) {
    return FAIL;
  }

  // compute FluxLine
  if (llf_mhd(FluxLine, priml, primr, prim, ActiveSize)==FAIL) {
    return FAIL;
  }

  if (NSpecies > 0) {
    plm_species(prim, 9, species, FluxLine[iD], ActiveSize);
    for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
	FluxLine[field][i] = FluxLine[iD][i]*species[field-NEQ_MHD][i];
      }
    }
  }

  if (NColor > 0) {
    plm_color(prim, 9, colors, FluxLine[iD], ActiveSize);
    for (int field = NEQ_MHD+NSpecies; field < NEQ_MHD+NSpecies+NColor; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
	FluxLine[field][i] = FluxLine[iD][i]*colors[field-NEQ_MHD-NSpecies][i];
      }
    }
  }
  
  return SUCCESS;
}
