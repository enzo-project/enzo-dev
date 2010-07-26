/***********************************************************************
/
/  HLL-PLM SOLVER
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
#include "ReconstructionRoutines.h"

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq);
int plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize);
int plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize);
int hll(float **FluxLine, float **priml, float **primr, int ActiveSize);

int HLL_PLM(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik)
{

  // compute priml and primr
  if (plm(prim, priml, primr, ActiveSize, 5) == FAIL) {
    return FAIL;
  }

  // compute FluxLine
  if (hll(FluxLine, priml, primr, ActiveSize)==FAIL) {
    return FAIL;
  }

  if (NSpecies > 0) {
    plm_species(prim, 5, species, FluxLine[iD], ActiveSize);
    for (int field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
	FluxLine[field][i] = FluxLine[iD][i]*species[field-NEQ_HYDRO][i];
      }
    }
  }

  if (NColor > 0) {
    plm_color(prim, 5, colors, FluxLine[iD], ActiveSize);
    for (int field = NEQ_HYDRO+NSpecies; field < NEQ_HYDRO+NSpecies+NColor; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
	FluxLine[field][i] = FluxLine[iD][i]*colors[field-NEQ_HYDRO-NSpecies][i];
      }
    }
  }
  
  return SUCCESS;
}
