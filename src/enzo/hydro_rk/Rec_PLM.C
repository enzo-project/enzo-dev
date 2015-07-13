/***********************************************************************
/
/  PLM RECONSTRUCTION
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

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ReconstructionRoutines.h"


inline void plm_l(float &vm1, float &v, float &vp1, float &vl_plm)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

inline void plm_r(float &vm1, float &v, float &vp1, float &vr_plm)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vr_plm = v - 0.5*dv;
}

inline void plm_point(float &vm1, float &v, float &vp1, float &vl_plm)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq)
{
  int iprim;

  for (int i = 0; i < ActiveSize+1; i++) {
    iprim = i + NumberOfGhostZones - 1;
    for (int field = 0; field < Neq; field++) {
      plm_point(prim[field][iprim-1], prim[field][iprim  ], prim[field][iprim+1],
		priml[field][i]);
      plm_point(prim[field][iprim+2], prim[field][iprim+1], prim[field][iprim],
		primr[field][i]);
    }

    priml[0][i] = max(priml[0][i], SmallRho);
    //priml[1][i] = max(priml[1][i], SmallEint);
    primr[0][i] = max(primr[0][i], SmallRho);
    //primr[1][i] = max(primr[1][i], SmallEint);
  }

  return SUCCESS;
}

int plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize)
  /* is : starting index of species field in prim */
{

  int iprim;
  for (int n = 0; n < ActiveSize+1; n++) {
    iprim = n + NumberOfGhostZones - 1;    
    if (flux0[n] >= 0) {
      for (int field = 0; field < NSpecies; field++) {
	plm_l(prim[field+is][iprim-1], prim[field+is][iprim], prim[field+is][iprim+1],
	      species[field][n]);
      }
    } else {
      for (int field = 0; field < NSpecies; field++) {
	plm_r(prim[field+is][iprim], prim[field+is][iprim+1], prim[field+is][iprim+2],
	      species[field][n]);
      }
    }

    /* renormalize species field */
    if (NoMultiSpeciesButColors != TRUE) {
      float sum = 0;
      for (int field = 0; field < NSpecies; field++) {
        sum += species[field][n];
      }
      for (int field = 0; field < NSpecies; field++) {
        species[field][n] /= sum;
      }
    }
  }

  return SUCCESS;
  
}

int plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize)
  /* is : starting index of *species* field in prim */
{

  int iprim;
  for (int n = 0; n < ActiveSize+1; n++) {
    iprim = n + NumberOfGhostZones - 1;        
    if (flux0[n] >= 0) {
      for (int field = is+NSpecies; field < is+NSpecies+NColor; field++) {
	plm_l(prim[field][iprim-1], prim[field][iprim], prim[field][iprim+1],
	      color[field-is-NSpecies][n]);
      }
    } else {
      for (int field = is+NSpecies; field < is+NSpecies+NColor; field++) {
	plm_r(prim[field][iprim], prim[field][iprim+1], prim[field][iprim+2],
	      color[field-is-NSpecies][n]);
      }
    }

  }

  return SUCCESS;
  
}
