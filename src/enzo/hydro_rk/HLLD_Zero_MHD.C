/***********************************************************************
/
/  HLLD WITH ZERO ORDER RECONSTRUCTION
/
/  written by: J. S. Oishi
/  date:       Apr, 2011
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

int hlld_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize);
int HLLD_Zero_MHD(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik)
{

  // compute priml and primr                                                                                
  int iprim;
  for (int field = 0; field < NEQ_MHD; field++) {
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i + NumberOfGhostZones - 1;
      priml[field][i] = prim[field][iprim];
      primr[field][i] = prim[field][iprim+1];
    }
  }

  // compute FluxLine
  if (hlld_mhd(FluxLine, priml, primr, prim, ActiveSize)==FAIL) {
    return FAIL;
  }

  float sum;
  if (NSpecies > 0) {
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i+NumberOfGhostZones-1;
      for (int field = 0; field < NSpecies; field++) {
        if (FluxLine[iD][i] >= 0) {
          species[field][i] = prim[field+5][iprim  ];
        } else {
          species[field][i] = prim[field+5][iprim+1];
        }
      }
      sum = 0;
      for (int field = 0; field < NSpecies; field++) {
        sum += species[field][i];
      }
      for (int field = 0; field < NSpecies; field++) {
        species[field][i] /= sum;
      }
      for (int field = 0; field < NSpecies; field++) {
        FluxLine[field+NEQ_MHD][i] = FluxLine[iD][i]*species[field][i];
      }
    }
  }

  if (NColor > 0) {
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i+NumberOfGhostZones-1;
      for (int field = 0; field < NColor; field++) {
        if (FluxLine[iD][i] >= 0) {
          colors[field][i] = prim[field+5+NSpecies][iprim  ];
        } else {
          colors[field][i] = prim[field+5+NSpecies][iprim+1];
        }
      }
      for (int field = 0; field < NColor; field++) {
        FluxLine[field+NEQ_MHD+NSpecies][i] = FluxLine[iD][i]*colors[field][i];
      }
    }
  }
  
  return SUCCESS;
}
