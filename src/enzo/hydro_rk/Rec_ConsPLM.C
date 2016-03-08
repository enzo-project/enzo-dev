/***********************************************************************
/
/  PLM RECONSTRUCTION
/
/  written by: Tom Abel
/  date:       October, 2009
/  modified1:
/
/ Conservative PLM Reconstruction
/ K Waagan
/ Journal of Computational Physics 228 (2009) 8609–8626 
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

inline void dv_plm_point(float &vm1, float &v, float &vp1, float &dv_plm)
{
  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv_plm = minmod(dv_l, dv_r, dv_m);

}

int cons_plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq, char direc)
{
  int iprim, iv;
  int idual = (DualEnergyFormalism) ? 1 : 0;
  float *dv[NEQ_MHD-idual];
  for (int field = 0; field < NEQ_MHD-idual; field++) 
    dv[field] = new float[ActiveSize+1];

  for (int i = 0; i < ActiveSize+1; i++) {
    iprim = i + NumberOfGhostZones - 1;
    for (int field = 0; field < Neq; field++) 
      dv_plm_point(prim[field][iprim-1], prim[field][iprim  ], prim[field][iprim+1],
		dv[field][i]);
  }

  // limiting certain slopes 3.14 
    for (int i = 0; i < ActiveSize+1; i++) {
      iprim = i + NumberOfGhostZones - 1;
      // density
      dv[0][i] = sign(dv[0][i])*min(fabs(dv[0][i]), prim[0][iprim]/2);
      // pressure
      dv[1][i] = sign(dv[1][i])*min(fabs(dv[1][i]), prim[1][iprim]/(1.+Gamma));
      // velocities
      for (int field = 2; field < 5; field++) 
	dv[field][i] = sign(dv[field][i])*min(fabs(dv[field][i]), 2./prim[9][iprim]) ;
    }


  /*    For U- reconstruction we would need more
     What's here now, Kaagan calles it p-reconstruction.  */

  // conservative reconstruction (equ.: 2.6-2.7)
  for (int i = 0; i < ActiveSize+1; i++) {
    iprim = i + NumberOfGhostZones - 1;
    // density & pressure
    
    for (int field = 0; field < 2; field++) { 
      priml[field][i] = prim[field][iprim] + 0.5*dv[field][i];
      primr[field][i] = prim[field][iprim] - 0.5*dv[field][i+1];
      if (field == 0) printf("field i :priml, prim, primr: %"ISYM" %"ISYM" : %"GSYM" %"GSYM" %"GSYM" \n", field, i, 
	     primr[field][i],prim[field][iprim], priml[field][i]);
    }
    // velocities
    for (int field = 2; field < 5; field++) {
      priml[field][i] = prim[field][iprim] + primr[0][i]/2/prim[0][iprim]*dv[field][i];
      primr[field][i] = prim[field][iprim] - priml[0][i]/2/prim[0][iprim]*dv[field][i+1];
    }
    // Bfields and phi
    for (int field = 5; field < Neq; field++) {
      priml[field][i] = prim[field][iprim] + 0.5*dv[field][i];
      primr[field][i] = prim[field][iprim] - 0.5*dv[field][i+1];
    }
  }

  return SUCCESS;

}
