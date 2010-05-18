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
/  Initial Guess Computation Routine
/
/  written by: Daniel Reynolds
/  date:       August, 2007
/  modified1:  
/
/  PURPOSE: Computes an initial guess to the time-evolved solution, 
/           according to the method specified by the initial_guess
/           input parameter:
/                0 -> use previous time step
/                1 -> full fwd Euler (local sources)
/                2 -> full fwd Euler (all rhs)
/                3 -> partial fwd Euler (local sources)
/                4 -> partial fwd Euler (all rhs)
/                5 -> use local analytical initial guess for reaction terms
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"

 
 

int gFLDProblem::InitialGuess(EnzoVector *uvec)
{
  // set some pointers
  float *tmp3_E, *Er, *Er0;

  // set initial guess depending on user choice
  switch (initial_guess) {


  case 0:  // previous time step

    uvec->copy(U0);
    break;


  case 1:  // full forward Euler step, local rhs sources only

    if (this->LocRHS(tmp3, told, U0) == FAIL) 
      ENZO_FAIL("InitialGuess Error: LocRHS failure");
    uvec->linearsum(1.0, U0, dt, tmp3);
    break;


  case 2:  // full forward Euler step, all rhs terms

    uvec->linearsum(1.0, U0, dt, rhs0);
    break;


  case 3:  // partial forward Euler step, local rhs sources only

    if (this->LocRHS(tmp3, told, U0) == FAIL) 
      ENZO_FAIL("InitialGuess Error: LocRHS failure");
    uvec->linearsum(1.0, U0, 0.1*dt, tmp3);
    break; 


  case 4:  // partial forward Euler step, all rhs terms

    uvec->linearsum(1.0, U0, 0.1*dt, rhs0);
    break;


  case 5:  // analytical initial guess for full reaction network

    //   Note: extsrc is available thanks to the last call to ComputeRHS.
    uvec->copy(U0);
    if (this->AnalyticInitGuess(uvec,dt) == FAIL) 
      ENZO_FAIL("InitialGuess Error: AnalyticInitGuess failure");
    break;


  default:  // illegal choice

    ENZO_FAIL("InitialGuess Error: illegal initial_guess choice");

  }

  return SUCCESS;

}

#endif
