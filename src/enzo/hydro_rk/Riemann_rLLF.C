/***********************************************************************
/
/  LLF RIEMANN SOLVER
/
/  written by: Tom Abel
/  date:       Feb, 2010
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

int rllf(float **FluxLine, float **priml, float **primr, int ActiveSize)
{
  float Ul[NEQ_HYDRO], Ur[NEQ_HYDRO], Fl[NEQ_HYDRO], Fr[NEQ_HYDRO];
  float etot, eintl, eintr, h, dpdrho, dpde, W, W2, ap, am, cs_l, cs_r, v_yz, v_yz2, v2,
    vx, vx2, vy, vz, rho, p, lm_l, lp_l, lm_r, lp_r, v, a0, pl, hl,rhol; 
  float Zero = 0.0;

  for (int n = 0; n < ActiveSize+1; n++) {
      // First, compute Fl and Ul
      rho   = priml[0][n];
      eintl = priml[1][n];
      vx    = priml[2][n];
      vy    = priml[3][n];
      vz    = priml[4][n];    
      
      // relative velocities only in the first pass
      
      EOS(p, rho, eintl, h, cs_l, dpdrho, dpde, EOSType, 2);
      if (EOSType > 0) {
	p = priml[1][n];
	pl = p;
	hl = h;
	cs_l = sqrt(p/rho);
      }

      Ul[iD   ] = rho;
      Ul[iS1  ] = 0.;
      Ul[iS2  ] = 0.;
      Ul[iS3  ] = 0.;
      Ul[iEtot] = 0.;
      if (DualEnergyFormalism) {
	Ul[iEint] = 0.;
      }
      
      Fl[iD   ] = 0.;
      Fl[iS1  ] = 0. + p  ; // Ul[iS1] * vx + p;
      Fl[iS2  ] = 0.; // Ul[iS2] * vx;
      Fl[iS3  ] = 0.; // Ul[iS3] * vx;
      Fl[iEtot] = 0.; // rho * (0.5*v2 + h) *vx;
      if (DualEnergyFormalism) {
	Fl[iEint] = 0.;
      }
      
      lp_l =   0;
      lm_l =   0;
      
      // Then, Fr and Ur
      rho   = primr[0][n];
      eintr = primr[1][n];
      vx    = primr[2][n];
      vy    = primr[3][n];
      vz    = primr[4][n];
      
      // relative velocity
      //      rho -=  priml[0][n];
      vx -=   priml[2][n];
      vy -=   priml[3][n];
      vz -=   priml[4][n];

      v2 = vx*vx + vy*vy + vz*vz;
      etot = eintr + 0.5*v2;
      EOS(p, rho, eintr, h, cs_r, dpdrho, dpde, EOSType, 2);
      if (EOSType > 0) {
	p = primr[1][n];
	cs_r = sqrt(p/rho);
      }
      
      rho -= rhol;
      Ur[iD   ] = rho;
      Ur[iS1  ] = rho * vx;
      Ur[iS2  ] = rho * vy;
      Ur[iS3  ] = rho * vz;
      Ur[iEtot] = rho * etot;
      if (DualEnergyFormalism) {
	Ur[iEint] = rho * eintr;
      }
      
      Fr[iD   ] = rho * vx;
      Fr[iS1  ] = Ur[iS1] * vx + p - pl;
      Fr[iS2  ] = Ur[iS2] * vx;
      Fr[iS3  ] = Ur[iS3] * vx;
      Fr[iEtot] = rho * (0.5*v2 + h - hl) *vx;
      if (DualEnergyFormalism) {
	Fr[iEint] = Ur[iEint] * vx;
      }
      
  

    ap = Max(0, lp_l, lp_r);
    am = Max(0, -lm_l, -lm_r);
    
    a0 = max(ap, am) + priml[2][n];

    for (int field = 0; field < NEQ_HYDRO; field++) {
      FluxLine[field][n] = 0.5*(Fr[field] - a0 * Ur[field]);
    }


    // now the actual left flux
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
    Ul[iEtot] = rho * etot;
    if (DualEnergyFormalism) {
      Ul[iEint] = rho * eintl;
    }

    Fl[iD   ] = rho * vx;
    Fl[iS1  ] = Ul[iS1] * vx ;
    Fl[iS2  ] = Ul[iS2] * vx;
    Fl[iS3  ] = Ul[iS3] * vx;
    Fl[iEtot] = rho * (0.5*v2) *vx;
    if (DualEnergyFormalism) {
      Fl[iEint] = Ul[iEint] * vx;
    }

    
    for (int field = 0; field < NEQ_HYDRO; field++) {
      FluxLine[field][n] += Fl[field];
    }


  }
  
  return SUCCESS;
}
