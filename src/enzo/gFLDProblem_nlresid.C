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
/  Gray Flux-Limited Diffusion Implicit Problem nonlinear residual 
/  function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Nonlinear residual function that defines the coupled,
/           implicit-time, radiation diffusion/chemistry/fluid energy 
/           system.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



int gFLDProblem::nlresid(EnzoVector *fu, EnzoVector *u)
{
#ifdef USE_JBPERF
    JBPERF_START("gfldproblem_nlresid");
#endif
  
//   if (debug)
//     printf("Entering gFLDProblem::nlresid routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) {
    fprintf(stderr,"nlresid error: gFLDProblem not yet prepared\n");
    return FAIL;
  }

//   // have u communicate neighbor information
//   if (u->exchange() == FAIL) {
//     fprintf(stderr,"nlresid error: EnzoVector::exchange failure\n");
//     return FAIL;
//   }
  // have u initiate communication of neighbor information
  if (u->exchange_start() == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::exchange_start failure\n");
    return FAIL;
  }

  // initialize residual to zero
  if (fu->constant(0.0) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::constant failure\n");
    return FAIL;
  }
  
  // have u finish communication of neighbor information
  if (u->exchange_end() == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::exchange_end failure\n");
    return FAIL;
  }

  // update time-dependent boundary conditions on state u, if needed
  if (this->UpdateBoundary(u,tnew,1) == FAIL) {
    fprintf(stderr,"nlresid error: UpdateBoundary failure\n");
    return FAIL;
  }

  // enforce boundary conditions on state u
  if (this->EnforceBoundary(u,0) == FAIL) {
    fprintf(stderr,"nlresid error: EnforceBoundary failure\n");
    return FAIL;
  }

  // compute rhs at current state for updated time
  if (this->ComputeRHS(rhs, tnew, u) == FAIL) {
    fprintf(stderr,"nlresid error: ComputeRHS failure\n");
    return FAIL;
  }

  // combine together via theta-scheme
  //   fu =  u-u0
  if (fu->linearsum(1.0,u,-1.0,U0) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::linearsum failure\n");
    return FAIL;
  }

  //   fu = (u-u0) - dt*(1-theta)*rhs0
  if (fu->axpy(-dt*(1.0-theta),rhs0)  == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::axpy failure\n");
    return FAIL;
  }

  //   fu = (u-u0) - dt*(1-theta)*rhs0 - dt*theta*rhs
  if (fu->axpy(-dt*theta,rhs) == FAIL) {
    fprintf(stderr,"nlresid error: EnzoVector::axpy failure\n");
    return FAIL;
  }

  
  // if using the analytical chemistry solver, call it now 
  // (overwrites fu components with analytical residual)
  int i;
  if (AnalyticChem == 1) {
    for (i=1; i<Nchem+2; i++)  fu->scale_component(i, 0.0);
    if (this->AnalyticResid(U0, u, fu, dt) == FAIL) {
      fprintf(stderr,"nlresid error: AnalyticalResid failure\n");
      return FAIL;
    }
  }


  //   Enforce boundary conditions on fu (for Dirichlet faces)
  if (this->EnforceBoundary(fu,1) == FAIL) {
    fprintf(stderr,"nlresid error: EnforceBoundary failure\n");
    return FAIL;    
  } 


//   // output component residuals
//   float rmsvals[Nchem+2];
//   for (i=0; i<Nchem+2; i++) 
//     rmsvals[i] = fu->infnorm_component(i);
//   if (debug) {
//     printf("    ");
//     for (i=0; i<Nchem+2; i++) 
//       printf("    f(%"ISYM") = %.2e",i,rmsvals[i]);
//     printf("\n");
//   }


#ifdef USE_JBPERF
    JBPERF_STOP("gfldproblem_nlresid");
#endif
  // return success
  return SUCCESS;

}
#endif
