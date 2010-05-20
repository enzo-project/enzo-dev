/***********************************************************************
/
/  HLL RIEMANN SOLVER
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

int hll(float **FluxLine, float **priml, float **primr, int ActiveSize)
{
  float Ul[NEQ_HYDRO], Ur[NEQ_HYDRO], Fl[NEQ_HYDRO], Fr[NEQ_HYDRO];
  float etot, eintl, eintr, h, dpdrho, dpde, W, W2, ap, am, cs_l, cs_r, v_yz, v_yz2, v2,
    vx, vx2, vy, vz, rho, p, lm_l, lp_l, lm_r, lp_r, v; 
  float Zero = 0.0;

  for (int n = 0; n < ActiveSize+1; n++) {
    // First, compute Fl and Ul
    rho   = priml[0][n];
    eintl = priml[1][n];
    vx    = priml[2][n];
    vy    = priml[3][n];
    vz    = priml[4][n];    

    v2 = vx*vx + vy*vy + vz*vz;
    etot = eintl + 0.5*v2;
    EOS(p, rho, eintl, h, cs_l, dpdrho, dpde, EOSType, 2);
    if (EOSType > 0) {
      p = priml[1][n];
      cs_l = sqrt(p/rho);
    }

    Ul[iD   ] = rho;
    Ul[iS1  ] = rho * vx;
    Ul[iS2  ] = rho * vy;
    Ul[iS3  ] = rho * vz;
    Ul[iEtot ] = rho * etot;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho * eintl;
    }

    Fl[iD   ] = rho * vx;
    Fl[iS1  ] = Ul[iS1] * vx + p;
    Fl[iS2  ] = Ul[iS2] * vx;
    Fl[iS3  ] = Ul[iS3] * vx;
    Fl[iEtot ] = rho * (0.5*v2 + h) *vx;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx;
    }

    lp_l = vx + cs_l;
    lm_l = vx - cs_l;

    // Then, Fr and Ur
    rho   = primr[0][n];
    eintr = primr[1][n];
    vx    = primr[2][n];
    vy    = primr[3][n];
    vz    = primr[4][n];

    v2 = vx*vx + vy*vy + vz*vz;
    etot = eintr + 0.5*v2;
    EOS(p, rho, eintr, h, cs_r, dpdrho, dpde, EOSType, 2);
    if (EOSType > 0) {
      p = primr[1][n];
      cs_r = sqrt(p/rho);
    }

    Ur[iD   ] = rho;
    Ur[iS1  ] = rho * vx;
    Ur[iS2  ] = rho * vy;
    Ur[iS3  ] = rho * vz;
    Ur[iEtot] = rho * etot;
    if (DualEnergyFormalism) {
      Ur[iEint] = rho * eintr;
    }

    Fr[iD   ] = rho * vx;
    Fr[iS1  ] = Ur[iS1] * vx + p;
    Fr[iS2  ] = Ur[iS2] * vx;
    Fr[iS3  ] = Ur[iS3] * vx;
    Fr[iEtot] = rho * (0.5*v2 + h) *vx;
    if (DualEnergyFormalism) {
      Fr[iEint] = Ur[iEint] * vx;
    }


    lp_r = vx + cs_r;
    lm_r = vx - cs_r;

    ap = Max(Zero, lp_l, lp_r);
    am = Max(Zero, -lm_l, -lm_r);

    for (int field = 0; field < NEQ_HYDRO; field++) {
      FluxLine[field][n] = (ap*Fl[field]+am*Fr[field]-ap*am*(Ur[field]-Ul[field]))/(ap+am);
    }

    /*if (fabs(FluxLine[iEint][n])>1e6) {
      printf("F[ieint] NaN at n=%"ISYM": fl[iS1]=%lf, fr[iS1]=%lf, Ul[iS1]=%lf, Ur[iS1]=%lf, ap=%lf, am = %lf\n",
	     n, Fl[iEint], Fr[iEint], Ul[iEint], Ur[iEint], ap, am);
      printf("lp_l=%lf, lm_l=%lf, lp_r=%lf, lm_r=%lf, cs_l=%lf, cs_r=%lf\n",
	     lp_l, lm_l, lp_r, lm_r, cs_l, cs_r);
      printf("priml: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     priml[0][n], priml[1][n], priml[2][n], priml[3][n], priml[4][n]);
      printf("primr: rho = %"GSYM", eint = %"GSYM", vx = %"GSYM", vy = %"GSYM", vz = %"GSYM"\n",
	     primr[0][n], primr[1][n], primr[2][n], primr[3][n], primr[4][n]);
      return FAIL;
      }*/
  }


  return SUCCESS;
}
