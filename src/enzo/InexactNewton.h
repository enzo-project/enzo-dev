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
/  Inexact Newton Solver Class
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: This class is designed to handle setup and solution of 
/           nonlinear systems using an inexact Newton method.
/
/           Implements a damped inexact Newton method:  Given an initial 
/           guess x0, while nonlinear residual is too large, do
/                  (i) solve J(x)*s = -f(x) inexactly
/                 (ii) update x += l*s, for 0 < l <= 1
/           where J(x) = d/dx[f(x)] is the Jacobian matrix.
/
/           The measure of nonlinear convergence is given via NtolNorm, 
/           where different values denote different measures of 
/           nonlinear convergence:
/                  0  ->  ||f(x)||_rms < Ntol
/                  1  ->  ||f(x)||_inf < Ntol
/                  2  ->  ||f(x)||_rms / ||f(x0)||_rms < Ntol
/                  3  ->  ||f(x)||_inf / ||f(x0)||_inf  < Ntol
/                  4  ->  ||f(x)||_rms / ||x||_rms < Ntol
/                  5  ->  ||f(x)||_inf / ||x||_inf < Ntol
/           where ||.||_rms denotes the root-mean-squared norm (2-norm 
/           divided by number of values) and ||.||_inf denotes the 
/           infinity norm (or maximum absolute value).  Choices 1 and 2 
/           require raw residual values below Ntol; choices 2 and 3 
/           require relative residual decrease below Ntol; choices 4 
/           and 5 weight the residual by the magnitude of the solution.
/           Control over the test type and Ntol are given through the
/           SetNewtonNorm() and SetNewtonTolerance() routines.
/
/           The inexactness parameter in step (i) is given through 
/           the InexactNewtonForce() routine.  The linesearch step size
/           l is found in the LinesearchStepSize() routine.  The linear 
/           solver used in step (i) is provided by the problem-specific 
/           lsolve() routine.
/
************************************************************************/

#ifndef INEXACT_NEWTON_SOLVER_DEFINED__
#define INEXACT_NEWTON_SOLVER_DEFINED__

#include "preincludes.h"
/* #include <stdlib.h> */
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "EnzoVector.h"
#include "NonlinearProblemABC.h"


class InexactNewtonSolver
{

 private:
  
  // temporary vectors for inexact Newton method
  EnzoVector *tmp1;
  EnzoVector *dx;
  EnzoVector *fvec;
  
  // solver parameters
  int    DampNewt;  // use damped Newton (linesearch) [!=0] or not [=0]
  int    Nmax_it;   // maximum number of Newton iterations
  int    INChoice;  // inexact Newton method
  int    NtolNorm;  // choice of norm for testing nonlinear convergence:
  float  Ntol;      // nonlinear solver tolerance
  float  IN1;       // inexact Newton parameter
  float  IN2;       // inexact Newton parameter
  float  LSmin;     // minimum line-search step size

  // solver diagnostics and internal constants 
  int    Niters;    // number of Newton iterations performed
  float  INforce;   // inexact Newton force
  float  fnorm;     // nonlinear residual norm
  float  fnorm0;    // temporary nonlinear residual norm
  float  LSsize;    // line-search step size


  // internal utility routines:
  //   computes the inexact Newton force
  int InexactNewtonForce();
  //   computes the linesearch step size
  int LinesearchStepSize(NonlinearProblemABC *prob, 
			 EnzoVector *x, EnzoVector *dx);


 public:


  // Constructor
  InexactNewtonSolver(EnzoVector *x);

  // Destructor
  ~InexactNewtonSolver();

  // Performs nonlinear solve: prob is Implicit Problem object, 
  //   x is vector object with initial guess (in) and solution (out)
  int Solve(NonlinearProblemABC *prob, EnzoVector *x);
  
  // Newton solver information routines:
  //   number of nonlinear iterations
  int GetNonlinearIterations();
  //   final nonlinear residual value
  float GetNonlinearResidual();
  //   final linesearch step size
  float GetLinesearchStepsize();

  // Newton solver set-parameter routine
  //   Maximum number of Newton iterations, must be >= 1
  int SetMaxIters(int maxit);
  //   Choice of Inexact Newton forcing method.  See 
  //   InexactNewton_InexactNewtonForce for options.
  int SetInexactNewton(int choice, float p1, float p2);
  //   Choice of nonlinear convergence measure, as described above.
  int SetNewtonNorm(int NormChoice);
  //   Nonlinear tolerance value, as described above.
  int SetNewtonTolerance(float tol);
  //   Minimum value for line-search step size, must be > 0.
  int SetMinLinesearch(float min);
  //   Damped Newton
  int SetDampedNewton(int choice);

};


#endif
