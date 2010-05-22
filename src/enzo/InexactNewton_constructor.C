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
/  Inexact Newton Solver constructor
/
/  written by: Daniel R. Reynolds
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:    Creates locally-owned temporary vectors, sets default 
/              solver parameters.
/
************************************************************************/
#include "performance.h"
#include "InexactNewton.h"


InexactNewtonSolver::InexactNewtonSolver(EnzoVector *x)
{
//   if (debug)
//     printf("Entering InexactNewtonSolver::constructor routine\n");

  // clone x to create tmp1, dx and fvec
  tmp1 = x->clone();
  dx   = x->clone();
  fvec = x->clone();
  
  // set default solver parameters
  DampNewt = 1;         // use the linesearch-enabled solver
  Nmax_it  = 20;        // maximum Newton iterations
  NtolNorm = 0;         // nonlinear tolerance norm (0->rms norm)
  Ntol     = 1.0e-7;    // nonlinear tolerance
  INChoice = 0;         // use constant inexactness forcing
  IN1      = 1.0e-1;    // inexactness forcing constant
  IN2      = 0.9;       // unused for now
  LSmin    = 1.0e-9;    // minimum line-search step size
};
