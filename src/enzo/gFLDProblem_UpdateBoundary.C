/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class 
/  UpdateBoundary routine
/
/  written by: Daniel Reynolds
/  date:       June, 2007
/  modified1:  
/
/  PURPOSE: Updates boundary conditions if necessary.  Depending on 
/           'flag' this routine will perform one of two BC-related 
/           actions:
/
/              flag=0: updates time-dependent boundary conditions 
/              for a given time step.  This will be called once per 
/              time step, during the 'Setup' phase.
/
/              flag!=0: updates solution-dependent boundary conditions.  
/              As these boundary conditions may change nonlinearly, due
/              to changes in the unknowns (gas energy, chemistry, etc.), 
/              this will be called at the beginning of the nonlinear 
/              residual calculation routine.  Hence, any solution-
/              dependent boundary conditions will only be treated using 
/              a simple fixed-point iteration, and not within the Newton 
/              system (since we do not include derivatives with respect 
/              to these BCs within the Newton solve).
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



int gFLDProblem::UpdateBoundary(EnzoVector *u, float time, int flag)
{
//   if (debug)
//     printf("Entering gFLDProblem::UpdateBoundary routine\n");

  // check that the gFLDProblem has been prepared
  if (!prepared) 
    ENZO_FAIL("UpdateBoundary ERROR: gFLDProblem unprepared");
  
  // get information about the vector u, and check against BC dims
  int i, j, k, l, m, idx, idxbc;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  if (udims[0] != LocDims[0]) {
    fprintf(stderr,"p%"ISYM" UpdateBC: mismatched x0 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[0],LocDims[0]);
    ENZO_FAIL("Error in gFLDProblem_UpdateBoundary");
  }
  if (udims[1] != LocDims[1]) {
    fprintf(stderr,"p%"ISYM" UpdateBC: mismatched x1 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[1],LocDims[1]);
    ENZO_FAIL("Error in gFLDProblem_UpdateBoundary");
  }
  if (udims[2] != LocDims[2]) {
    fprintf(stderr,"p%"ISYM" UpdateBC: mismatched x2 dims %"ISYM"!=%"ISYM"\n",
	    MyProcessorNumber,udims[2],LocDims[2]);
    ENZO_FAIL("Error in gFLDProblem_UpdateBoundary");
  }
  if (udims[3] != (2+Nchem)) {
    fprintf(stderr,"p%"ISYM" UpdateBC: mismatched nspecies %"ISYM"!=3\n",
	    MyProcessorNumber,udims[3]);
    ENZO_FAIL("Error in gFLDProblem_UpdateBoundary");
  }

  // set some shortcuts for the EnzoVector dimensions
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];


  // perform boundary condition updates based on Model

  //   Marshak test problems: update time-dependent BCs
  if ( Model >= 20 && Model <= 29 ) {

    // first update all nonlinearly-evolving BCs
    // (if put here, we perform a fixed-point 
    //  iteration on the nonlinearity in the BC fluxes)

    // next update remaining time-dependent BCs 
    if (flag == 0) {

      // set some relevant constants
      float bval, kap, tau, c, pi, so_eps;
      idx = ((ugh[2][0])*x1len+ugh[1][0])*x0len;  // first interior cell index
      kap = OpacityE[idx];                        // opacity (constant)
      so_eps = MarshakParms[0];                   // 'epsilon' for Su & Olson
      c = 2.99792458e10;                          // speed of light (cm/s)
      pi = 4.0*atan(1.0);                         // no explanation needed
      tau = so_eps * c * kap * time;              // non-dimensional S&O 'time'

      // compute the time-dependent boundary value:
      // if tau is small, use asymptotic solution [Su & Olson, eq (51)], 
      if (tau < 1e-5) {
	bval = sqrt(3.0*tau/pi/so_eps) - 3.0*tau/4/so_eps 
	     + tau/2.0/so_eps*sqrt(tau/3.0/pi/so_eps);
	//printf("gFLDProblem_UpdateBoundary: bval = %g\n",bval);
      }
      // otherwise use M-level Romberg numerical integration of 
      // analytical solution [Su & Olson, eq (36)-(38) after x=0 manipulations]
      else {
	float integral1, integral2;
	int M=19;
	float R1[M+1][M+1], R2[M+1][M+1];
	float a=1.0e-15, b=1.0-a, h, x, g1, g2, f1, f2;
	for (k=0; k<M+1; k++)
	  for (l=0; l<M+1; l++) {
	    R1[k][l] = 0.0;
	    R2[k][l] = 0.0;
	  }
	
	//   compute first Romberg entry using midpoint rule
	h = b-a;
	x = a;
	g1 = x*sqrt(so_eps + 1.0/(1.0-x*x));
	f1 = exp(-tau*x*x)*(2.0*g1)/(x*(3.0+4.0*g1*g1));
	R1[0][0] = 0.5*h*f1;
	g2 = sqrt((1.0-x)*(so_eps+1.0/x));
	f2 = exp(-tau/so_eps/x)*(2.0*g2)/(x*(1.0+so_eps*x)*(3.0+4.0*g2*g2));
	R2[0][0] = 0.5*h*f2;
	
	x = b;
	g1 = x*sqrt(so_eps + 1.0/(1.0-x*x));
	f1 = exp(-tau*x*x)*(2.0*g1)/(x*(3.0+4.0*g1*g1));
	R1[0][0] += 0.5*h*f1;
	g2 = sqrt((1.0-x)*(so_eps+1.0/x));
	f2 = exp(-tau/so_eps/x)*(2.0*g2)/(x*(1.0+so_eps*x)*(3.0+4.0*g2*g2));
	R2[0][0] += 0.5*h*f2;
	
	//   iterate over remaining Romberg rows
	for (k=1; k<=M; k++) {
	  // adjust mesh spacing for this level
	  h = 0.5*h;
	  
	  // compute first entry of next Romberg row
	  R1[k][0] = 0.0;
	  R2[k][0] = 0.0;
	  for (l=1; l<=POW(2,k-1); l++) {
	    x = a - h + 2.0*h*l;
	    g1 = x*sqrt(so_eps + 1.0/(1.0-x*x));
	    f1 = exp(-tau*x*x)*(2.0*g1)/(x*(3.0+4.0*g1*g1));
	    R1[k][0] += f1;
	    g2 = sqrt((1.0-x)*(so_eps+1.0/x));
	    f2 = exp(-tau/so_eps/x)*(2.0*g2)/(x*(1.0+so_eps*x)*(3.0+4.0*g2*g2));
	    R2[k][0] += f2;
	  } // end l loop
	  R1[k][0] = h*R1[k][0] + 0.5*R1[k-1][0];
	  R2[k][0] = h*R2[k][0] + 0.5*R2[k-1][0];
	  
	  // compute remaining entries of Romberg row
	  for (m=1; m<=k; m++) {
	    R1[k][m] = (R1[k][m-1]-R1[k-1][m-1])/(POW(4.0,m)-1.0) + R1[k][m-1];
	    R2[k][m] = (R2[k][m-1]-R2[k-1][m-1])/(POW(4.0,m)-1.0) + R2[k][m-1];
	  } // end m loop
	} // end k loop
	
	//   set integrals to last entries in Romberg tables
	integral1 = R1[M][M];
	integral2 = R2[M][M];
	
	// set boundary value according to analytical solution
	bval = 1.0 - (sqrt(3.0)/pi) * (2.0*integral1 + exp(-tau)*integral2);
      }

      // fill x0 left Dirichlet boundary
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idxbc = k*LocDims[1] + j;
	  EBdryVals[0][0][idxbc] = bval/ErUnits;
	}
    }

  } // end if (20 <= Model < 30)

  // return success
  return SUCCESS;
}
#endif
