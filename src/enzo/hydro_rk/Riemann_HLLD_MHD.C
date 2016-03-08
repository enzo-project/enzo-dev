/***********************************************************************
/
/  HLLD RIEMANN SOLVER
/
/  written by: J. S. Oishi
/  date:       1 April 2011
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

int hlld_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize)
{
  float Ul[NEQ_MHD], Ur[NEQ_MHD], Fl[NEQ_MHD], Fr[NEQ_MHD], Us[NEQ_MHD], Uss[NEQ_MHD];
  float etot_l,etot_r, eint_l, eint_r, h, dpdrho, dpde, rho_l, rho_r, vx_l, vy_l, vz_l, vx_r, vy_r, vz_r, Bx_l, Bx_r,Bx, By_l, Bz_l, By_r, Bz_r, Phi_l, Phi_r, v2, B2, Bv_l, Bv_r, p_l, p_r, cs_l, cs_r, pt_l, pt_r;
  float rho_ls, rho_rs, vy_ls, vy_rs, vz_ls, vz_rs, vv_ls, vv_rs, By_ls, By_rs, Bz_ls, Bz_rs, Bv_ls, Bv_rs, bb_ls, bb_rs, eint_ls, eint_rs, etot_ls, etot_rs, pt_s;
  float vy_ss, vz_ss, By_ss, Bz_ss, Bv_ss, eint_lss, eint_rss, etot_lss, etot_rss, rho_savg;
  float Zero = 0.0;
  float S_l, S_r, S_ls, S_rs, S_M; // wave speeds
  float cf_l, cf_r, sam, sap; // fast speeds

  for (int n = 0; n < ActiveSize+1; n++) {
    // First, compute Fl and Ul
    rho_l  = priml[0][n];
    eint_l = priml[1][n];
    vx_l   = priml[2][n];
    vy_l   = priml[3][n];
    vz_l   = priml[4][n];    
    Bx_l   = priml[5][n];
    By_l   = priml[6][n];
    Bz_l   = priml[7][n];
    Phi_l   = priml[8][n];
    B2 = Bx_l * Bx_l + By_l * By_l + Bz_l * Bz_l;
    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;

    v2 = vx_l * vx_l + vy_l * vy_l + vz_l * vz_l;
    etot_l = rho_l * (eint_l + 0.5 * v2) + 0.5 * B2;
    EOS(p_l, rho_l, eint_l, h, cs_l, dpdrho, dpde, EOSType, 2);
    pt_l = p_l + 0.5 * B2;
    cf_l = sqrt((Gamma * p_l + B2 + sqrt((Gamma * p_l + B2) * (Gamma * p_l + B2) - 4. * Gamma * p_l * Bx_l * Bx_l))/(2. * rho_l));

    Ul[iD   ] = rho_l;
    Ul[iS1  ] = rho_l * vx_l;
    Ul[iS2  ] = rho_l * vy_l;
    Ul[iS3  ] = rho_l * vz_l;
    Ul[iEtot] = etot_l;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho_l * eint_l;
    }
    Ul[iBx ] = Bx_l;
    Ul[iBy ] = By_l;
    Ul[iBz ] = Bz_l;
    Ul[iPhi] = Phi_l;

    Fl[iD   ] = rho_l * vx_l;
    Fl[iS1  ] = Ul[iS1] * vx_l + pt_l - Bx_l * Bx_l;
    Fl[iS2  ] = Ul[iS2] * vx_l - Bx_l * By_l;
    Fl[iS3  ] = Ul[iS3] * vx_l - Bx_l * Bz_l;
    Fl[iEtot] = (etot_l + pt_l) * vx_l - Bx_l * Bv_l;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx_l;
    }
    Fl[iBx] = 0.0;
    Fl[iBy] = vx_l*By_l - vy_l*Bx_l;
    Fl[iBz] = -vz_l*Bx_l + vx_l*Bz_l;

    //compute Ur and Fr
    rho_r   = primr[0][n];
    eint_r  = primr[1][n];
    vx_r    = primr[2][n];
    vy_r    = primr[3][n];
    vz_r    = primr[4][n];
    Bx_r    = primr[5][n];
    By_r    = primr[6][n];
    Bz_r    = primr[7][n];
    Phi_r   = primr[8][n];
    B2 = Bx_r * Bx_r + By_r * By_r + Bz_r * Bz_r;
    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;

    v2 = vx_r * vx_r + vy_r * vy_r + vz_r * vz_r;
    etot_r = rho_r * (eint_r + 0.5 * v2) + 0.5 * B2;
    EOS(p_r, rho_r, eint_r, h, cs_r, dpdrho, dpde, EOSType, 2);
    pt_r = p_r + 0.5 * B2;
    cf_r = sqrt((Gamma * p_r + B2 + sqrt((Gamma * p_r + B2) * (Gamma * p_r + B2) - 4. * Gamma * p_r * Bx_r * Bx_r))/(2. * rho_r));

    Ur[iD   ] = rho_r;
    Ur[iS1  ] = rho_r * vx_r;
    Ur[iS2  ] = rho_r * vy_r;
    Ur[iS3  ] = rho_r * vz_r;
    Ur[iEtot] = etot_r;
    if (DualEnergyFormalism) {
      Ur[iEint] = rho_r * eint_r;
    }
    Ur[iBx ] = Bx_r;
    Ur[iBy ] = By_r;
    Ur[iBz ] = Bz_r;
    Ur[iPhi] = Phi_r;

    Fr[iD   ] = rho_r * vx_r;
    Fr[iS1  ] = Ur[iS1] * vx_r + pt_r - Bx_r * Bx_r;
    Fr[iS2  ] = Ur[iS2] * vx_r - Bx_r * By_r;
    Fr[iS3  ] = Ur[iS3] * vx_r - Bx_r * Bz_r;
    Fr[iEtot] = (etot_r + pt_r) * vx_r - Bx_r * Bv_r;
    if (DualEnergyFormalism) {
      Fr[iEint] = Ur[iEint] * vx_l;
    }
    Fr[iBx ] = 0.0;
    Fr[iBy ] = vx_r * By_r - vy_r * Bx_r;
    Fr[iBz ] = -vz_r * Bx_r + vx_r * Bz_r;

    //
    //wave speeds
    //

    Bx = 0.5*(Bx_l + Bx_r);
    // first, outermost wave speeds
    // simplest choice from Miyoshi & Kusano (2005)
    S_l = min(vx_l, vx_r) - max(cf_l, cf_r);
    S_r = max(vx_l, vx_r) + max(cf_l, cf_r);

    if (S_l > 0) {
      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fl[field];
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);

      continue;
    } 
    if (S_r < 0) {
      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fr[field];
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);

      continue;
    } 

    // next, the middle (contact) wave
    S_M = ((S_r - vx_r)*rho_r*vx_r - (S_l - vx_l)*rho_l*vx_l - pt_r + pt_l)/((S_r - vx_r)*rho_r - (S_l - vx_l)*rho_l);

    // finally, the intermediate (Alfven) waves

    rho_ls = rho_l * (S_l - vx_l)/(S_l - S_M);
    rho_rs = rho_r * (S_r - vx_r)/(S_r - S_M);

    S_ls = S_M - fabs(Bx)/sqrt(rho_ls);
    S_rs = S_M + fabs(Bx)/sqrt(rho_rs);

    pt_s = ((S_r -  vx_r) * rho_r*pt_l - (S_l - vx_l) * rho_l * pt_r + rho_l*rho_r*(S_r - vx_r)*(S_l - vx_l)*(vx_r - vx_l))/((S_r - vx_r)*rho_r - (S_l - vx_l)*rho_l);

    sam = vx_l - cf_l;
    sap = vx_l + cf_l;
      
    if ((fabs(S_M - vx_l) <= BFLOAT_EPSILON) and 
        (fabs(By_l) <= BFLOAT_EPSILON) and 
        (fabs(Bz_l) <= BFLOAT_EPSILON) and 
        (Bx*Bx >= Gamma * p_l) and
        ((fabs(S_l - sam) <= BFLOAT_EPSILON) or (fabs(S_l - sap) <= BFLOAT_EPSILON)) ) {
      vy_ls = vy_l;
      vz_ls = vz_l;
      By_ls = By_l;
      Bz_ls = Bz_l;
    } else {
      vv_ls = (S_M - vx_l)/(rho_l*(S_l - vx_l)*(S_l - S_M) - Bx*Bx);
      bb_ls = (rho_l*(S_l - vx_l)*(S_l - vx_l) - Bx*Bx)/(rho_l*(S_l - vx_l)*(S_l - S_M) - Bx*Bx);
      vy_ls = vy_l - Bx * By_l * vv_ls;
      By_ls = By_l * bb_ls;
      vz_ls = vz_l - Bx * Bz_l * vv_ls;
      Bz_ls = Bz_l * bb_ls;
    }

    sam = vx_r - cf_r;
    sap = vx_r + cf_r;
      
    if ((fabs(S_M - vx_r) <= BFLOAT_EPSILON) and 
        (fabs(By_r) <= BFLOAT_EPSILON) and 
        (fabs(Bz_r) <= BFLOAT_EPSILON) and 
        (Bx*Bx >= Gamma * p_r) and
        ((fabs(S_r - sam) <= BFLOAT_EPSILON) or (fabs(S_r - sap) <= BFLOAT_EPSILON)) ) {
      vy_rs = vy_r;
      vz_rs = vz_r;
      By_rs = By_r;
      Bz_rs = Bz_r;
    } else {
      vv_rs = (S_M - vx_r)/(rho_r*(S_r - vx_r)*(S_r - S_M) - Bx*Bx);
      bb_rs = (rho_r*(S_r - vx_r)*(S_r - vx_r) - Bx*Bx)/(rho_r*(S_r - vx_r)*(S_r - S_M) - Bx*Bx);
      vy_rs = vy_r - Bx * By_r * vv_rs;
      vz_rs = vz_r - Bx * Bz_r * vv_rs;
      By_rs = By_r * bb_rs;
      Bz_rs = Bz_r * bb_rs;
    }
    Bv_ls = S_M * Bx + vy_ls * By_ls + vz_ls * Bz_ls;
    Bv_rs = S_M * Bx + vy_rs * By_rs + vz_rs * Bz_rs;

    etot_ls = ((S_l - vx_l)*etot_l - pt_l*vx_l + pt_s * S_M + Bx*(Bv_l - Bv_ls))/(S_l - S_M);
    etot_rs = ((S_r - vx_r)*etot_r - pt_r*vx_r + pt_s * S_M + Bx*(Bv_r - Bv_rs))/(S_r - S_M);
    
    // compute the fluxes based on the wave speeds
    if (S_l <= 0 && S_ls >= 0) {
      // USE F_ls
      Us[iD   ] = rho_ls;
      Us[iS1  ] = rho_ls * S_M;
      Us[iS2  ] = rho_ls * vy_ls;
      Us[iS3  ] = rho_ls * vz_ls;
      Us[iEtot] = etot_ls;
      if (DualEnergyFormalism)
        Us[iEint] = rho_ls * eint_ls;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_ls;
      Us[iBz  ] = Bz_ls;
      Us[iPhi ] = Phi_l;

      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fl[field] + S_l*(Us[field] - Ul[field]);
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);

      continue;
    } 
    if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      Us[iD   ] = rho_rs;
      Us[iS1  ] = rho_rs * S_M;
      Us[iS2  ] = rho_rs * vy_rs;
      Us[iS3  ] = rho_rs * vz_rs;
      Us[iEtot] = etot_rs;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_rs;
      Us[iBz  ] = Bz_rs;
      Us[iPhi ] = Phi_r;

      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fr[field] + S_r*(Us[field] - Ur[field]);
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);

      continue;
    } 

    //do U** stuff
    rho_savg = sqrt(rho_ls) + sqrt(rho_rs);
    vy_ss = (sqrt(rho_ls) * vy_ls + sqrt(rho_rs) * vy_rs + (By_rs - By_ls) * sign(Bx))/rho_savg;
    vz_ss = (sqrt(rho_ls) * vz_ls + sqrt(rho_rs) * vz_rs + (Bz_rs - Bz_ls) * sign(Bx))/rho_savg;
    By_ss = (sqrt(rho_ls) * By_rs + sqrt(rho_rs) * By_ls + sqrt(rho_ls * rho_rs) * (vy_rs - vy_ls) * sign(Bx))/rho_savg;
    Bz_ss = (sqrt(rho_ls) * Bz_rs + sqrt(rho_rs) * Bz_ls + sqrt(rho_ls * rho_rs) * (vz_rs - vz_ls) * sign(Bx))/rho_savg;
    Bv_ss = S_M * Bx + vy_ss * By_ss + vz_ss * Bz_ss;
    etot_lss = etot_ls - sqrt(rho_ls) * (Bv_ls - Bv_ss) * sign(Bx);
    etot_rss = etot_rs + sqrt(rho_rs) * (Bv_rs - Bv_ss) * sign(Bx);

    if (S_ls <= 0 && S_M >= 0) {
      // USE F_lss
      Us[iD   ] = rho_ls;
      Us[iS1  ] = rho_ls * S_M;
      Us[iS2  ] = rho_ls * vy_ls;
      Us[iS3  ] = rho_ls * vz_ls;
      Us[iEtot] = etot_ls;
      if (DualEnergyFormalism)
        Us[iEint] = rho_ls * eint_ls;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_ls;
      Us[iBz  ] = Bz_ls;
      Us[iPhi ] = Phi_l;

      Uss[iD   ] = rho_ls;
      Uss[iS1  ] = rho_ls * S_M;
      Uss[iS2  ] = rho_ls * vy_ss;
      Uss[iS3  ] = rho_ls * vz_ss;
      Uss[iEtot] = etot_lss;
      if (DualEnergyFormalism)
        Uss[iEint] = rho_ls * eint_lss;
      Uss[iBx  ] = Bx;
      Uss[iBy  ] = By_ss;
      Uss[iBz  ] = Bz_ss;
      Uss[iPhi ] = Phi_l;
      
      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fl[field] + S_ls*Uss[field] - (S_ls - S_l)*Us[field] - S_l*Ul[field];
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);
      continue;
    } 
    if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss
      Us[iD   ] = rho_rs;
      Us[iS1  ] = rho_rs * S_M;
      Us[iS2  ] = rho_rs * vy_rs;
      Us[iS3  ] = rho_rs * vz_rs;
      Us[iEtot] = etot_rs;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_rs;
      Us[iBz  ] = Bz_rs;
      Us[iPhi ] = Phi_r;

      Uss[iD   ] = rho_rs;
      Uss[iS1  ] = rho_rs * S_M;
      Uss[iS2  ] = rho_rs * vy_ss;
      Uss[iS3  ] = rho_rs * vz_ss;
      Uss[iEtot] = etot_rss;
      if (DualEnergyFormalism)
        Uss[iEint] = rho_rs * eint_rss;
      Uss[iBx  ] = Bx;
      Uss[iBy  ] = By_ss;
      Uss[iBz  ] = Bz_ss;
      Uss[iPhi ] = Phi_r;
      
      for (int field = 0; field < NEQ_MHD - 1; field++) {
	FluxLine[field][n] = Fr[field] + S_rs*Uss[field] - (S_rs - S_r)*Us[field] - S_r*Ur[field];
      }
      FluxLine[iBx][n]  = Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
      FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
      FluxLine[iPhi][n] *= (C_h*C_h);
      continue;
    }
  }
  return SUCCESS;
}
