/***********************************************************************
/
/  PPM RECONSTRUCTION
/
/  written by: Peng Wang 
/  date:       May, 2007
/  modified1:  Tom Abel adopted form relativistic version of the code 
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

void ppm_quartic(float **prim, float **p0l, float **p0r, int ActiveSize, int Neq)
  // Input: prim[Neq][ActiveSize+6]
  // Output: p0l,r[Neq+1][ActiveSize+3] 
  // note: the first element of p0l and the last element of p0r will not be used
  // in the future.
{

  float dap1, da, ap1mam1, amam1, ap1ma, ap2ma, ap2map1, a, ap2, ap1, am1;
  int iprim;

  for (int field = 0; field < Neq; field++) {
    for (int i = 0; i < ActiveSize+2; i++) {
      iprim = i + 1;

      ap2 = prim[field][iprim+2];
      ap1 = prim[field][iprim+1];
      a = prim[field][iprim];
      am1 = prim[field][iprim-1];
      // compute da_j+1
      ap2ma = ap2 - a;
      ap1ma = ap1 - a;
      ap2map1 = ap2 - ap1;

      if (ap2map1*ap1ma > 0) {
	dap1 = sign(ap2ma)*Min(0.5*fabs(ap2ma), 2.0*fabs(ap1ma), 2.0*fabs(ap2map1));
      } else {
	dap1 = 0.0;
      }

      // compute da_j
      ap1mam1 = ap1 - am1;
      amam1 = a - am1;
      ap1ma = ap1 - a;

      if (ap1ma*amam1 > 0) {
	da = sign(ap1mam1)*Min(0.5*fabs(ap1mam1), 2.0*fabs(amam1), 2.0*fabs(ap1ma));
      } else {
	da = 0.0;
      }

      p0l[field][i] = a + 0.5*ap1ma - dap1/6.0 + da/6.0;
      p0r[field][i] = p0l[field][i];
    }
  }
}

void ppm_dissipation(float **prim, float **p0l, float **p0r, int ActiveSize, int Neq)
  // prim[Neq][ActiveSize+6]
  // p0l,r[Neq][ActiveSize+3]
  // Reference: Colella & Woodward (1984), Marti & Muller (1996)  
{

  float K0 = 1.0, eta1 = 5.0, eta2 = 0.05, eps1 = 0.1, eps2 = 1.0, w1 = 0.52, w2 = 10.0;
  float dap1, dam1, amam2, am1mam2, amam1, ap2ma, ap1ma, ap2map1, ap1mam1, a, ap2, ap1, am1, am2,
    d2a, d2ap1, d2am1, eta0, eta, adr, adl;
  float ar, al;
  float f0, f0p1, f0m1, f, p_p3, p_m3; 
  float rho_p1, rho_m1, p_p1, p_m1, v_p1, v_m1, p_p2, p_m2, p;

  int iprim;
  int ipres=1;
  for (int i = 0; i < ActiveSize+1; i++) {
    iprim = i + 2;
    rho_p1 = prim[iden][iprim+1];
    rho_m1 = prim[iden][iprim-1];
    p = prim[ipres][iprim];
    p_p1 = prim[ipres][iprim+1];
    p_p2 = prim[ipres][iprim+2];
    p_m1 = prim[ipres][iprim-1];
    p_m2 = prim[ipres][iprim-2];
    v_p1 = prim[ivx][iprim+1];
    v_m1 = prim[ivx][iprim-1];

    for (int field = 0; field < Neq; field++) {
      ap2 = prim[field][iprim+2];
      ap1 = prim[field][iprim+1];
      a = prim[field][iprim];
      am1 = prim[field][iprim-1];
      am2 = prim[field][iprim-2];
     
      if (1) {
	// 1) Contact steeping
	if (Gamma*K0*fabs(rho_p1-rho_m1)/min(rho_p1,rho_m1) >= 
	    fabs(p_p1-p_m1)/min(p_p1, p_m1)) {
	  
	  // compute da_j+1
	  ap2ma = ap2 - a;
	  ap1ma = ap1 - a;
	  ap2map1 = ap2 - ap1;
	  
	  if (ap2map1*ap1ma > 0) {
	    dap1 = sign(ap2ma)*Min(0.5*fabs(ap2ma), 2.0*fabs(ap1ma), 2.0*fabs(ap2map1));
	  } else {
	    dap1 = 0.0;
	  }
	  
	  // compute da_j-1
	  amam2 = a - am2;
	  am1mam2 = am1 - am2;
	  amam1 = a - am1;
	  
	  if (amam1*am1mam2 > 0) {
	    dam1 = sign(amam2)*Min(0.5*fabs(amam2), 2.0*fabs(am1mam2), 2.0*fabs(amam1));
	  } else {
	    dam1 = 0.0;
	  }
	  
	  adl = ap1 - 0.5*dap1;
	  adr = am1 + 0.5*dam1;
	  
	  // compute eta_j
	  d2a = ap1 - 2.0*a + am1;
	  d2ap1 = ap2 - 2.0*ap1 + a;
	  d2am1 = a - 2.0*am1 + am2;
	  ap1mam1 = ap1 - am1;
	  if (-d2ap1*d2am1 > 0 && fabs(ap1mam1) - eps1*min(fabs(ap1),fabs(am1)) > 0) {
	    eta0 = -(d2a-d2am1)/(6.0*ap1mam1);
	  } else {
	    eta0 = 0.0;
	  }
	  
	  eta = max(0.0, min(eta1*(eta0-eta2), 1.0));
	  
	  p0l[field][i+1] = p0l[field][i+1]*(1.0-eta) + adl*eta;
	  p0r[field][i] = p0r[field][i]*(1.0-eta) + adr*eta;
	}
	
	// 2) Shock flattening
	if (i >= 1 && i <= ActiveSize) {
	  if (fabs(p_p1-p_m1)/min(p_p1,p_m1) > eps2 && v_m1 > v_p1) {
	    f0 = min(1.0, max(0.0, ((p_p1-p_m1)/(p_p2-p_m2) - w1)*w2));
	    if (p_p1 - p_m1 > 0) {
	      p_p3 = prim[ipres][iprim+3];
	      f0p1 = min(1.0, max(0.0, ((p_p2-p)/(p_p3-p_m1) - w1)*w2));
	      f = max(f0, f0p1);
	    } else {
	      p_m3 = prim[ipres][iprim-3];
	      f0m1 = min(1.0, max(0.0, ((p-p_m2)/(p_p1-p_m3) - w1)*w2));
	      f = max(f0, f0m1);
	    }
	    
	    p0l[field][i+1] = a*f + p0l[field][i+1]*(1.0-f);
	    p0r[field][i] = a*f + p0r[field][i]*(1.0-f);
	  }
	}
      } // turn on/off steepening and shock flattening

      // 3) Monotonization:
      ar = p0l[field][i+1];
      al = p0r[field][i];
      if ((ar-a)*(a-al) <= 0) {
	p0l[field][i+1] = a;
	p0r[field][i] = a;
      }
      if ((ar-al)*(a-0.5*(ar+al)) > (ar-al)*(ar-al)/6.0) {
	p0r[field][i] = 3.0*a - 2.0*ar;
      }
      if (-(ar-al)*(a-0.5*(ar+al)) > (ar-al)*(ar-al)/6.0) {
	p0l[field][i+1] = 3.0*a - 2.0*al;
      }
    }   
  }
}

int ppm(float **prim, float **priml, float **primr, int ActiveSize, int Neq)
{
  ppm_quartic(prim, priml, primr, ActiveSize, Neq);
  ppm_dissipation(prim, priml, primr, ActiveSize, Neq);
  return SUCCESS;
}
