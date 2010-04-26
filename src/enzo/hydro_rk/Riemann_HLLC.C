/***********************************************************************
/
/  HLLC RIEMANN SOLVER
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
#include "EOS.h"  

int hllc(float **FluxLine, float **priml, float **primr, int ActiveSize)
{
  float Ul[NEQ_HYDRO], Ur[NEQ_HYDRO], Fl[NEQ_HYDRO], Fr[NEQ_HYDRO], Uc[NEQ_HYDRO], Fc[NEQ_HYDRO];
  float etot, eint_l, eint_r, h, dpdrho, dpde, W, W2, ap, am, cs_l, cs_r, v_yz, v_yz2, v2,
    vx_l, vx_r, vx2, vy_l, vy_r, vz_r, vz_l, rho_r, rho_l, p_l, p_r, lm_l, lp_l, lm_r, lp_r, v; 
  float eint_c, rho_c, p_c, vx_c, vy_c, vz_c, cs_c;
  float Zero = 0.0, lam_l, lam_r, lam_c;

  for (int n = 0; n < ActiveSize+1; n++) {
    // First, compute Fl and Ul
    rho_l  = priml[0][n];
    eint_l = priml[1][n];
    vx_l   = priml[2][n];
    vy_l   = priml[3][n];
    vz_l   = priml[4][n];    

    v2 = vx_l*vx_l + vy_l*vy_l + vz_l*vz_l;
    etot = eint_l + 0.5*v2;
    EOS(p_l, rho_l, eint_l, h, cs_l, dpdrho, dpde, EOSType, 2);

    Ul[iD   ] = rho_l;
    Ul[iS1  ] = rho_l * vx_l;
    Ul[iS2  ] = rho_l * vy_l;
    Ul[iS3  ] = rho_l * vz_l;
    Ul[iEtot] = rho_l * etot;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho_l * eint_l;
    }

    Fl[iD   ] = rho_l * vx_l;
    Fl[iS1  ] = Ul[iS1] * vx_l + p_l;
    Fl[iS2  ] = Ul[iS2] * vx_l;
    Fl[iS3  ] = Ul[iS3] * vx_l;
    Fl[iEtot] = rho_l * (0.5*v2 + h) *vx_l;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx_l;
    }

    lp_l = vx_l + cs_l;
    lm_l = vx_l - cs_l;

    // Then, Fr and Ur
    rho_r  = primr[0][n];
    eint_r = primr[1][n];
    vx_r   = primr[2][n];
    vy_r   = primr[3][n];
    vz_r   = primr[4][n];

    v2 = vx_r*vx_r + vy_r*vy_r + vz_r*vz_r;
    etot = eint_r + 0.5*v2;
    EOS(p_r, rho_r, eint_r, h, cs_r, dpdrho, dpde, EOSType, 2);

    Ur[iD   ] = rho_r;
    Ur[iS1  ] = rho_r * vx_r;
    Ur[iS2  ] = rho_r * vy_r;
    Ur[iS3  ] = rho_r * vz_r;
    Ur[iEtot] = rho_r * etot;
    if (DualEnergyFormalism) {
      Ur[iEint] = rho_r * eint_r;
    }

    Fr[iD   ] = rho_r * vx_r;
    Fr[iS1  ] = Ur[iS1] * vx_r + p_r;
    Fr[iS2  ] = Ur[iS2] * vx_r;
    Fr[iS3  ] = Ur[iS3] * vx_r;
    Fr[iEtot] = rho_r * (0.5*v2 + h) * vx_r;
    if (DualEnergyFormalism) {
      Fr[iEint] = Ur[iEint] * vx_r;
    }


    lp_r = vx_r + cs_r;
    lm_r = vx_r - cs_r;

    lam_l = min(lm_l, lm_r);
    lam_r = max(lp_l, lp_r);


    lam_c = (rho_r*vx_r*(lam_r-vx_r) - rho_l*vx_l*(lam_l-vx_l)+p_l-p_r)/
      (rho_r*(lam_r-vx_r) - rho_l*(lam_l-vx_l));

    if (lam_l >= 0) {
      for (int field = 0; field < NEQ_HYDRO; field++) {
	FluxLine[field][n] = Fl[field];
      }
    } else if (lam_r <= 0) {
      for (int field = 0; field < NEQ_HYDRO; field++) {
	FluxLine[field][n] = Fr[field];
      }
    } else if (lam_c >= 0) {
      // compute Fc_l
      rho_c = rho_l*(lam_l-vx_l)/(lam_l-lam_c);
      p_c = p_l + rho_l*(vx_l-lam_l)*(vx_l-lam_c);
      vx_c = lam_c;
      vy_c = vy_l*rho_l/rho_c*(lam_l-vx_l)/(lam_l-lam_c);
      vz_c = vz_l*rho_l/rho_c*(lam_l-vx_l)/(lam_l-lam_c);
      v2 = vx_c*vx_c + vy_c*vy_c + vz_c*vz_c;
      EOS(p_c, rho_c, eint_c, h, cs_c, dpdrho, dpde, EOSType, 1);
      Fc[iD  ] = rho_c * vx_c;
      Fc[iS1 ] = rho_c * vx_c * vx_c + p_c;
      Fc[iS2 ] = rho_c * vy_c * vx_c;
      Fc[iS3 ] = rho_c * vz_c * vx_c;
      Fc[iEtot] = rho_c * (0.5*v2 + h) * vx_c;
      if (DualEnergyFormalism) {
	Fc[iEint] = rho_c * eint_c * vx_c;
      }
      for (int field = 0; field < NEQ_HYDRO; field++) {
	FluxLine[field][n] = Fc[field];
      }
    } else {
      // compute Fc_r
      rho_c = rho_r*(lam_r-vx_r)/(lam_r-lam_c);
      p_c = p_r + rho_r*(vx_r-lam_r)*(vx_r-lam_c);
      vx_c = lam_c;
      vy_c = vy_r*rho_r/rho_c*(lam_r-vx_r)/(lam_r-lam_c);
      vz_c = vz_r*rho_r/rho_c*(lam_r-vx_r)/(lam_r-lam_c);
      v2 = vx_c*vx_c + vy_c*vy_c + vz_c*vz_c;
      EOS(p_c, rho_c, eint_c, h, cs_c, dpdrho, dpde, EOSType, 1);
      Fc[iD  ] = rho_c * vx_c;
      Fc[iS1 ] = rho_c * vx_c * vx_c + p_c;
      Fc[iS2 ] = rho_c * vy_c * vx_c;
      Fc[iS3 ] = rho_c * vz_c * vx_c;
      Fc[iEtot] = rho_c * (0.5*v2 + h) * vx_c;
      if (DualEnergyFormalism) {
	Fc[iEint] = rho_c * eint_c * vx_c;
      }
      for (int field = 0; field < NEQ_HYDRO; field++) {
	FluxLine[field][n] = Fc[field];
      }
    }
    

    if (isnan(FluxLine[iD][n])) {
      printf("F[iS1] NaN at n=%"ISYM": fl[iS1]=%lf, fr[iS1]=%lf, Ul[iS1]=%lf, Ur[iS1]=%lf, ap=%lf, am = %lf\n",
	     n, Fl[iD], Fr[iD], Ul[iD], Ur[iD], ap, am);
      printf("lp_l=%lf, lm_l=%lf, lp_r=%lf, lm_r=%lf, cs_l=%lf, cs_r=%lf\n",
	     lp_l, lm_l, lp_r, lm_r, cs_l, cs_r);
      printf("priml: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     priml[0][n], priml[1][n], priml[2][n], priml[3][n], priml[4][n]);
      printf("primr: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     primr[0][n], primr[1][n], primr[2][n], primr[3][n], primr[4][n]);
      return FAIL;
    }
  }


  /*for (int n = 0; n < ActiveSize; n++) {
    if (fabs(FluxLine[iEtot][n+1]-FluxLine[iEtot][n]) > 1e5) {
      printf("flux error\n");
      for (int n = 0; n < ActiveSize+1; n++) {
	printf("%"GSYM" ", FluxLine[iEtot][n]);
      }
      printf("\n");

      for (int n = 0; n < ActiveSize+1; n++) {
	rho  = priml[0][n];
	float eintl = priml[1][n];
	vx   = priml[2][n];
	vy   = priml[3][n];
	vz   = priml[4][n];    
	
	v2 = vx*vx + vy*vy + vz*vz;
	etot = eintl + 0.5*v2;
	EOS(p, rho, eintl, h, cs_l, dpdrho, dpde, 0, 2);
	printf("pl=%"GSYM", cs_l=%"GSYM" entl=%"GSYM", rho=%"GSYM" ",p, cs_l, eintl, rho);
	
	Ul[iD  ] = rho;
	Ul[iS1 ] = rho * vx;
	Ul[iS2 ] = rho * vy;
	Ul[iS3 ] = rho * vz;
	Ul[iEtot] = rho * etot;
	
	Fl[iD  ] = rho * vx;
	Fl[iS1 ] = Ul[iS1] * vx + p;
	Fl[iS2 ] = Ul[iS2] * vx;
	Fl[iS3 ] = Ul[iS3] * vx;
	Fl[iEtot] = rho * (0.5*v2 + h) *vx;
	
	lp_l = vx + cs_l;
	lm_l = vx - cs_l;
	
	// Then, Fr and Ur
	rho  = primr[0][n];
	float eintr = primr[1][n];
	vx   = primr[2][n];
	vy   = primr[3][n];
	vz   = primr[4][n];
	
	v2 = vx*vx + vy*vy + vz*vz;
	etot = eintr + 0.5*v2;
	EOS(p, rho, eintr, h, cs_r, dpdrho, dpde, 0, 2);
	printf("pr=%"GSYM",cs_r=%"GSYM" eintr=%"GSYM" rhor=%"GSYM" \n", p, cs_r, eintr, rho);
	
	Ur[iD  ] = rho;
	Ur[iS1 ] = rho * vx;
	Ur[iS2 ] = rho * vy;
	Ur[iS3 ] = rho * vz;
	Ur[iEtot] = rho * etot;
	
	Fr[iD  ] = rho * vx;
	Fr[iS1 ] = Ur[iS1] * vx + p;
	Fr[iS2 ] = Ur[iS2] * vx;
	Fr[iS3 ] = Ur[iS3] * vx;
	Fr[iEtot] = rho * (0.5*v2 + h) *vx;
    
	lp_r = vx + cs_r;
	lm_r = vx - cs_r;
    
	ap = Max(Zero, lp_l, lp_r);
	am = Max(Zero, -lm_l, -lm_r);
	printf("fr=%"GSYM",fl=%"GSYM",ap=%"GSYM",am=%"GSYM", ur=%"GSYM", ul=%"GSYM"\n",Fr[iEtot],Fl[iEtot],ap,am, Ur[iEtot], Ul[iEtot]);
	
	for (int field = 0; field < NEQ_HYDRO; field++) {
	  FluxLine[field][n] = (ap*Fl[field]+am*Fr[field]-ap*am*(Ur[field]-Ul[field]))/(ap+am);
	}
	
	
      }
   
      
      return FAIL;
    }
    }*/


  //return FAIL;
  return SUCCESS;
}
