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
//   if (debug)
//     printf("Entering gFLDProblem::nlresid routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) 
    ENZO_FAIL("nlresid error: gFLDProblem not yet prepared");

//   // have u communicate neighbor information
//   if (u->exchange() == FAIL) 
//      ENZO_FAIL("nlresid error: EnzoVector::exchange failure");
  // have u initiate communication of neighbor information
  if (u->exchange_start() == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::exchange_start failure");

  // initialize residual to zero
  if (fu->constant(0.0) == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::constant failure");
  
  // have u finish communication of neighbor information
  if (u->exchange_end() == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::exchange_end failure");

  // update time-dependent boundary conditions on state u, if needed
  if (this->UpdateBoundary(u,tnew,1) == FAIL) 
    ENZO_FAIL("nlresid error: UpdateBoundary failure");

  // enforce boundary conditions on state u
  if (this->EnforceBoundary(u,0) == FAIL) 
    ENZO_FAIL("nlresid error: EnforceBoundary failure");

  // compute rhs at current state for updated time
  if (this->ComputeRHS(rhs, tnew, u) == FAIL) 
    ENZO_FAIL("nlresid error: ComputeRHS failure");

  // combine together via theta-scheme
  //   fu =  u-u0
  if (fu->linearsum(1.0,u,-1.0,U0) == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::linearsum failure");

  //   fu = (u-u0) - dt*(1-theta)*rhs0
  if (fu->axpy(-dt*(1.0-theta),rhs0)  == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::axpy failure");

  //   fu = (u-u0) - dt*(1-theta)*rhs0 - dt*theta*rhs
  if (fu->axpy(-dt*theta,rhs) == FAIL) 
    ENZO_FAIL("nlresid error: EnzoVector::axpy failure");

  
  // if using the analytical chemistry solver, call it now 
  // (overwrites fu components with analytical residual)
  int i;
  if (AnalyticChem == 1) {
    for (i=1; i<Nchem+2; i++)  fu->scale_component(i, 0.0);
    if (this->AnalyticResid(U0, u, fu, dt) == FAIL) 
      ENZO_FAIL("nlresid error: AnalyticalResid failure");
  }


  //   Enforce boundary conditions on fu (for Dirichlet faces)
  if (this->EnforceBoundary(fu,1) == FAIL) 
    ENZO_FAIL("nlresid error: EnforceBoundary failure");


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

  // return success
  return SUCCESS;

}
#endif
