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

int hlld(float **FluxLine, float **priml, float **primr, int ActiveSize)
{
  float Ul[NEQ_MHD], Ur[NEQ_MHD], Fl[NEQ_MHD], Fr[NEQ_MHD], Fc[NEQ_MHD];
  float etot, eint_l, eint_r, h, dpdrho, dpde, rho_l, eint_l, vx_l, vy_l, vz_l, Bx, By, Bz, Phi;

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
    Phi_l  = priml[8][n];
    B2 = Bx*Bx + By*By + Bz*Bz;
    Bv_l = Bx*vx + By*vy + Bz*vz;

    v2 = vx*vx + vy*vy + vz*vz;
    etot = eint + 0.5*v2 + 0.5*B2/rho;
    EOS(p_l, rho_l, eint_l, h, cs_l, dpdrho, dpde, EOSType, 2);
    pt_l = p_l + 0.5*B2;
    Ul[iD   ] = rho_l;
    Ul[iS1  ] = rho_l * vx_l;
    Ul[iS2  ] = rho_l * vy_l;
    Ul[iS3  ] = rho_l * vz_l;
    Ul[iEtot] = rho_l * etot;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho_l * eint_l;
    }
    Ul[iBx ] = Bx;
    Ul[iBy ] = By;
    Ul[iBz ] = Bz;
    Ul[iPhi] = Phi;

    Fl[iD   ] = rho_l * vx_l;
    Fl[iS1  ] = Ul[iS1] * vx_l + p_l;
    Fl[iS2  ] = Ul[iS2] * vx_l;
    Fl[iS3  ] = Ul[iS3] * vx_l;
    Fl[iEtot] = rho_l * (0.5*v2 + h) *vx_l + B2*vx_l - Bx*Bv_l;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx_l;
    }
    Fl[iBx] = 0.0;
    Fl[iBy] = vx*By - vy*Bx;
    Fl[iBz] = -vz*Bx + vx*Bz;

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

    B2 = Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r;
    Bv_r = Bx_r*vx_r + By_r*vy_r + Bz_r*vz_r;

    v2 = vx_r*vx_r + vy_r*vy_r + vz_r*vz_r;
    etot_r = eint_r + 0.5*v2 + 0.5*B2/rho_r;
    EOS(p_r, rho_r, eint_r, h, cs_r, dpdrho, dpde, EOSType, 2);
    pt_r = p_r + 0.5*B2;
    // if (EOSType > 0) {
    //   p = primr[1][n];
    //   cs = sqrt(p/rho);
    // }
    //float pr=p;
    //cs2 = cs*cs;

    Ur[iD   ] = rho_r;
    Ur[iS1  ] = rho_r * vx_r;
    Ur[iS2  ] = rho_r * vy_r;
    Ur[iS3  ] = rho_r * vz_r;
    Ur[iEtot] = rho_r * etot_r;
    if (DualEnergyFormalism) {
      Ur[iEint] = rho_r * eint_r;
    }
    Ur[iBx ] = Bx;
    Ur[iBy ] = By_r;
    Ur[iBz ] = Bz_r;
    Ur[iPhi] = Phi_r;

    Fr[iD   ] = rho_r * vx_r;
    Fr[iS1  ] = Ur[iS1] * vx_r + p_r + 0.5*B2 - Bx*Bx;
    Fr[iS2  ] = Ur[iS2] * vx_r - Bx*By_r;
    Fr[iS3  ] = Ur[iS3] * vx_r - Bx*Bz_r;
    Fr[iEtot] = rho_r * (0.5*v2 + h) *vx_r + B2*vx_r - Bx*Bv_r;
    if (DualEnergyFormalism) {
      Fr[iEint] = Ur[iEint] * vx;
    }
    Fr[iBx ] = 0.0;
    Fr[iBy ] = vx*By_r - vy_r*Bx;
    Fr[iBz ] = -vz*Bx + vx_r*Bz_r;
    //Fr[iPhi] = C_h*C_h*Bx;

    //
    //wave speeds
    //

    // first, outermost wave speeds
    // These are not right! Fix!
    c_p = vx_r + cs_r;  
    c_m = vx_r - cs_r;
    S_l = Max(Zero, c_p, c_m);
    c_p = vx_l + cs_l;
    c_l = vx_l - cs_l;
    S_r = Max(Zero, -c_p, -c_m);

    if (S_l > 0) {
      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fl[field];
      }
      continue;
    } 
    if (S_r < 0) {
      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fr[field];
      }
      continue;
    } 

    // next, the middle (contact) wave
    S_M = ((S_r - vx_r)*rho_r*vx_r - (S_l - vx_l)*rho_l*vx_l - p_tr + p_tl)/((S_r * vx_r)*rho_r - (S_l - vx_l)*rho_l);

    // finally, the intermediate (Alfven) waves
    rho_ls = (S_l * vx_l - rho_l * vx_l)/(S_l = S_M);
    rho_rs = (S_r * vx_r - rho_r * vx_r)/(S_r = S_M);
    S_ls = S_M - abs(Bx)/sqrt(rho_ls);
    S_rs = S_M + abs(Bx)/sqrt(rho_rs);

    pt_s = ((S_r -  vx_r) * rho_r*pt_l - (S_l - vx_l) * rho_l * pt_r + rho_l*rho_r*(S_r - vx_r)*(S_l - vx_l)*(vx_r - vx_l))/((S_r - vx_r)*rho_r - (S_l - u_l)*rho_l);

    vv_ls = (S_M - vx_l)/(rho_l*(S_l - vx_l)*(S_l - S_M) - Bx*Bx);
    vv_rs = (S_M - vx_r)/(rho_r*(S_r - vx_r)*(S_r - S_M) - Bx*Bx);
    bb_lxs = (rho_l*(S_l - vx_l)*(S_l - vx_l) - B_x*B_x)/(rho_l*(S_l - vx_l)*(S_l - S_M) - B_x*B_x);
    bb_lys = (rho_r*(S_r - vx_r)*(S_r - vx_r) - B_x*B_x)/(rho_r*(S_r - vx_r)*(S_r - S_M) - B_x*B_x);

    vy_ls = vy_l - Bx * By_l * vv_ls;
    vy_rs = vy_r - Bx * By_r * vv_rs;
    By_ls = By_l*bb_lxs;
    By_rs = By_r*bb_lys;

    vz_ls = vz_l - Bx * Bz_l * vv_ls;
    vz_rs = vz_r - Bx * Bz_l * vv_rs;
    Bz_ls = Bz_l * bb_ls;
    Bz_rs = Bz_r * bb_rs;

    e_ls = ((S_l * vx_l)*e_l - pt_l*vx_l + pt_s * S_M + Bx*(Bv_l - Bv_ls))/(S_l - S_M);
    e_rs = ((S_r * vx_r)*e_r - pt_r*vx_r + pt_s * S_M + Bx*(Bv_r - Bv_rs))/(S_r - S_M);
    
    // compute the fluxes based on the wave speeds
    if (S_l <= 0 && S_ls >= 0) {
      // USE F_ls
      Us[iD   ] = rho_ls;
      Us[iS1  ] = rho_ls * vx_l;
      Us[iS2  ] = rho_ls * vy_l;
      Us[iS3  ] = rho_ls * vz_l;
      Us[iEtot] = rho_ls * e_ls;
      if (DualEnergyFormalism)
        Us[iEint] = rho_ls * eint_ls;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_ls;
      Us[iBz  ] = Bz_ls;
      Us[iPhi ] = Phi;

      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fl[i] + S_l*(Us[i] - Ul[i]);
      }
      continue;
    } 
    if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      Us[iD   ] = rho_rs;
      Us[iS1  ] = rho_rs * S_M;
      Us[iS2  ] = rho_rs * vy_r;
      Us[iS3  ] = rho_rs * vz_r;
      Us[iEtot] = rho_rs * e_rs;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_rs;
      Us[iBz  ] = Bz_rs;
      Us[iPhi ] = Phi;

      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fr[i] + S_r*(Us[i] - Ul[i]);
      }
      
      continue;
    } 

    //do U** stuff
    vy_ss = sqrt(rho_ls) * vy_ls + sqrt(rho_rs) * vy_rs + (By_rs - By_ls) * sign(Bx);
    vz_ss = sqrt(rho_ls) * vz_ls + sqrt(rho_rs) * vz_rs + (Bz_rs - Bz_ls) * sign(Bx);
    By_ss = sqrt(rho_ls) * By_rs + sqrt(rho_rs)  * By_ls + sqrt(rho_ls * rho_rs) * (vy_rs - vy_ls) * sign(Bx);
    Bz_ss = sqrt(rho_ls) * Bz_rs + sqrt(rho_rs) * Bz_ls + sqrt(rho_ls * rho_rs) * (vz_rs - vz_ls) * sign(Bx);
    Bv_ss = vx * Bx + vy_ss * By_ss + vz_ss * Bz_ss;
    etot_lss = etot_ls - sqrt(rho_ls) * (Bv_ls - Bv_ss) * sign(Bx);
    etot_rss = etot_rs - sqrt(rho_rs) * (Bv_rs - Bv_ss) * sign(Bx);

    if (S_ls <= 0 && S_M >= 0) {
      // USE F_lss
      Us[iD   ] = rho_ls;
      Us[iS1  ] = rho_ls * vx_l;
      Us[iS2  ] = rho_ls * vy_l;
      Us[iS3  ] = rho_ls * vz_l;
      Us[iEtot] = rho_ls * e_ls;
      if (DualEnergyFormalism)
        Us[iEint] = rho_ls * eint_ls;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_ls;
      Us[iBz  ] = Bz_ls;
      Us[iPhi ] = Phi;

      Uss[iD   ] = rho_ls;
      Uss[iS1  ] = rho_ls * S_M;
      Uss[iS2  ] = rho_ls * vy_ss;
      Uss[iS3  ] = rho_ls * vz_ss;
      Uss[iEtot] = rho_ls * etot_lss;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Uss[iBx  ] = Bx;
      Uss[iBy  ] = By_ss;
      Uss[iBz  ] = Bz_ss;
      Uss[iPhi ] = Phi;
      
      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fl[field] + S_ls*Uss[field] - (S_ls - S_l)*Us[field]*S_l*Ul[i];
      }
      
      continue;
    } 
    if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss

      // USE F_lss
      Us[iD   ] = rho_rs;
      Us[iS1  ] = rho_rs * vx_r;
      Us[iS2  ] = rho_rs * vy_r;
      Us[iS3  ] = rho_rs * vz_r;
      Us[iEtot] = rho_rs * e_rs;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Us[iBx  ] = Bx;
      Us[iBy  ] = By_rs;
      Us[iBz  ] = Bz_rs;
      Us[iPhi ] = Phi;

      Uss[iD   ] = rho_rs;
      Uss[iS1  ] = rho_rs * S_M;
      Uss[iS2  ] = rho_rs * vy_ss;
      Uss[iS3  ] = rho_rs * vz_ss;
      Uss[iEtot] = rho_rs * etot_rss;
      if (DualEnergyFormalism)
        Us[iEint] = rho_rs * eint_rs;
      Uss[iBx  ] = Bx;
      Uss[iBy  ] = By_ss;
      Uss[iBz  ] = Bz_ss;
      Uss[iPhi ] = Phi;
      
      for (int field = 0; field < NEQ_MHD; field++) {
	FluxLine[field][n] = Fr[field] + S_rs*Uss[field] - (S_rs - S_r)*Us[field]*S_r*Ur[i];
      }

      continue;
    }
    
  }
  return SUCCESS;
}
