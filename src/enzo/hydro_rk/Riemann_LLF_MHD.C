/***********************************************************************
/
/  MHD LLF RIEMANN SOLVER
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

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "../hydro_rk/ReconstructionRoutines.h"
#include "EOS.h"  

int llf_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize)
{
  float Ul[NEQ_MHD], Ur[NEQ_MHD], Fl[NEQ_MHD], Fr[NEQ_MHD];
  float etot, eint, h, dpdrho, dpde, W, W2, ap, am, cs, cs2, ca2, ca, cf, cf2, v_yz, v_yz2, v2,
    vx, vx2, vy, vz, rho, p, lm_l, lp_l, lm_r, lp_r, v, Bx, By, Bz, Phi, B2, Bv; 
  float Zero = 0.0;
  float temp1;

  for (int n = 0; n < ActiveSize+1; n++) {
    // First, compute Fl and Ul
    rho   = priml[0][n];
    eint  = priml[1][n];
    vx    = priml[2][n];
    vy    = priml[3][n];
    vz    = priml[4][n];
    Bx    = priml[5][n];
    By    = priml[6][n];
    Bz    = priml[7][n];
    Phi   = priml[8][n];
    
    B2 = Bx*Bx + By*By + Bz*Bz;
    Bv = Bx*vx + By*vy + Bz*vz;

    v2 = vx*vx + vy*vy + vz*vz;
    etot = eint + 0.5*v2 + 0.5*B2/rho;

    EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);    
    if (EOSType > 0) {
      p = priml[1][n];
      cs = sqrt(p/rho);
    }
    float pl = p;
    cs2 = cs*cs;

    Ul[iD   ] = rho;
    Ul[iS1  ] = rho * vx;
    Ul[iS2  ] = rho * vy;
    Ul[iS3  ] = rho * vz;
    Ul[iEtot] = rho * etot;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho * eint;
    }
    Ul[iBx ] = Bx;
    Ul[iBy ] = By;
    Ul[iBz ] = Bz;
    Ul[iPhi] = Phi;

    Fl[iD   ] = rho * vx;
    Fl[iS1  ] = Ul[iS1] * vx + p + 0.5*B2 - Bx*Bx;
    Fl[iS2  ] = Ul[iS2] * vx - Bx*By;
    Fl[iS3  ] = Ul[iS3] * vx - Bx*Bz;
    Fl[iEtot] = rho * (0.5*v2 + h) *vx + B2*vx - Bx*Bv;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx;
    }
    Fl[iBx] = 0.0;
    Fl[iBy] = vx*By - vy*Bx;
    Fl[iBz] = -vz*Bx + vx*Bz;
    //Fl[iPhi] = C_h*C_h*Bx;

    // largest and smallest eigenvectors
    ca2 = Bx*Bx/rho;
    temp1 = cs2 + B2/rho;
    cf2 = 0.5 * (temp1 + sqrt(fabs(temp1*temp1 - 4.0*cs2*ca2)));
    cf = sqrt(cf2);

    lp_l = vx + cf;
    lm_l = vx - cf;

    // Then, Fr and Ur
    rho   = primr[0][n];
    eint  = primr[1][n];
    vx    = primr[2][n];
    vy    = primr[3][n];
    vz    = primr[4][n];
    Bx    = primr[5][n];
    By    = primr[6][n];
    Bz    = primr[7][n];
    Phi   = primr[8][n];

    B2 = Bx*Bx + By*By + Bz*Bz;
    Bv = Bx*vx + By*vy + Bz*vz;

    v2 = vx*vx + vy*vy + vz*vz;
    etot = eint + 0.5*v2 + 0.5*B2/rho;
    EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
    if (EOSType > 0) {
      p = primr[1][n];
      cs = sqrt(p/rho);
    }
    float pr=p;
    cs2 = cs*cs;

    Ur[iD   ] = rho;
    Ur[iS1  ] = rho * vx;
    Ur[iS2  ] = rho * vy;
    Ur[iS3  ] = rho * vz;
    Ur[iEtot] = rho * etot;
    if (DualEnergyFormalism) {
      Ur[iEint] = rho * eint;
    }
    Ur[iBx ] = Bx;
    Ur[iBy ] = By;
    Ur[iBz ] = Bz;
    Ur[iPhi] = Phi;

    Fr[iD   ] = rho * vx;
    Fr[iS1  ] = Ur[iS1] * vx + p + 0.5*B2 - Bx*Bx;
    Fr[iS2  ] = Ur[iS2] * vx - Bx*By;
    Fr[iS3  ] = Ur[iS3] * vx - Bx*Bz;
    Fr[iEtot] = rho * (0.5*v2 + h) *vx + B2*vx - Bx*Bv;
    if (DualEnergyFormalism) {
      Fr[iEint] = Ur[iEint] * vx;
    }
    Fr[iBx ] = 0.0;
    Fr[iBy ] = vx*By - vy*Bx;
    Fr[iBz ] = -vz*Bx + vx*Bz;
    //Fr[iPhi] = C_h*C_h*Bx;

    // largest and smallest eigenvectors
    ca2 = Bx*Bx/rho;
    temp1 = cs2 + B2/rho;
    cf2 = 0.5 * (temp1 + sqrt(fabs(temp1*temp1 - 4.0*cs2*ca2)));
    cf = sqrt(cf2);

    lp_r = vx + cf;
    lm_r = vx - cf;

    ap = Max(Zero, lp_l, lp_r);
    am = Max(Zero, -lm_l, -lm_r);
    //ap = max(Zero, lp_r);
    //am = max(Zero, -lm_l);

    float a0 = max(ap, am);

    for (int field = 0; field < NEQ_MHD-1; field++) {
      FluxLine[field][n] = 0.5*(Fl[field]+Fr[field]-a0*(Ur[field]-Ul[field]));
    }


    FluxLine[iBx][n] += Ul[iPhi] + 0.5*(Ur[iPhi]-Ul[iPhi]) - 0.5*C_h*(Ur[iBx]-Ul[iBx]);
    FluxLine[iPhi][n] = Ul[iBx] + 0.5*(Ur[iBx]-Ul[iBx]) - 0.5/C_h*(Ur[iPhi]-Ul[iPhi]);
    FluxLine[iPhi][n] *= (C_h*C_h);

    /*if (isnan(FluxLine[iD][n])) {
      printf("F[iD] huge at n=%"ISYM": fl[iS1]=%lf, fr[iS1]=%lf, Ul[iS1]=%lf, Ur[iS1]=%lf, ap=%lf, am = %lf\n",
	     n, Fl[iD], Fr[iD], Ul[iD], Ur[iD], ap, am);
      printf("lp_l=%lf, lm_l=%lf, lp_r=%lf, lm_r=%lf, pl=%"GSYM", pr=%"GSYM", cfr=%"GSYM", temp2=%"GSYM", 4.0csca=%"GSYM"\n",
	     lp_l, lm_l, lp_r, lm_r, pl, pr, cf, temp1*temp1, 4.0*cs2*ca2);
      printf("cs2=%"GSYM", ca2=%"GSYM", B2/rho=%"GSYM"\n", cs2, ca2, B2/rho);
      printf("priml: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     priml[0][n], priml[1][n], priml[2][n], priml[3][n], priml[4][n]);
      printf("primr: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     primr[0][n], primr[1][n], primr[2][n], primr[3][n], primr[4][n]);
      return FAIL;
      }*/
  }

  /*for (int n = 0; n < ActiveSize+1; n++) {
    FluxLine[iBx][n] = 0.0;
    }*/

  return SUCCESS;
}
