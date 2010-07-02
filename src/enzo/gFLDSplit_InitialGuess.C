/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  Initial Guess Computation Routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Computes an initial guess to the time-evolved solution, 
/           according to the method specified by the initial_guess
/           input parameter:
/                 0 -> use previous time step
/              else -> use local analytical initial guess
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

 
 

int gFLDSplit::InitialGuess(EnzoVector *uvec)
{
  // set initial guess depending on user choice
  switch (initial_guess) {


  case 0:  // previous time step

    uvec->copy(U0);
    break;


  case 1:
  case 2:
  case 3:
  case 4:
  case 5:  // analytical initial guess for full reaction network
    //   Note: extsrc is available thanks to the last call to ComputeRHS.
    uvec->copy(U0);
    if (this->AnalyticInitGuess(uvec,dt) == FAIL) 
      ENZO_FAIL("InitialGuess Error: AnalyticInitGuess failure");
    break;


  default:  // illegal choice

    fprintf(stderr,"InitialGuess Error: illegal initial_guess choice = %"ISYM"\n",
	    initial_guess);
    ENZO_FAIL("Error in gFLDSplit_InitialGuess");

  }

  return SUCCESS;

}

#endif
