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
#include "InexactNewton.h"


// Newton solver information routines:

//   number of nonlinear iterations
int InexactNewtonSolver::GetNonlinearIterations() {
  return Niters;
}

//   final nonlinear residual value
float InexactNewtonSolver::GetNonlinearResidual() {
  return fnorm;
}

//   final linesearch step size
float InexactNewtonSolver::GetLinesearchStepsize() {
  return LSsize;
}



// Newton solver set-parameter routines

//   Maximum number of Newton iterations, must be >= 1
int InexactNewtonSolver::SetMaxIters(int maxit) {
  if (maxit > 0) { 
    Nmax_it = maxit; 
    return SUCCESS; 
  }
  else  
    return FAIL;
}

//   Choice of Inexact Newton forcing method.  See 
//   InexactNewton_InexactNewtonForce for options.
int InexactNewtonSolver::SetInexactNewton(int choice, float p1, float p2) {
  INChoice = choice;  
  IN1 = p1;  
  IN2 = p2;  
  return SUCCESS;
}

//   Choice of nonlinear convergence measure, as described above.
int InexactNewtonSolver::SetNewtonNorm(int NormChoice) {
  if ((NormChoice < 0) || (NormChoice > 5))  {
    NtolNorm = 0;
    fprintf(stderr," SetNewtonNorm warning: %"ISYM" out of bounds\n",NormChoice);
  }
  else  NtolNorm = NormChoice;
  return SUCCESS;
}

//   Nonlinear tolerance value, as described above.
int InexactNewtonSolver::SetNewtonTolerance(float tol) {
  if (tol > 0.0)  {
    Ntol = tol;  
    return SUCCESS; 
  }
  else  
    return FAIL;
}

//   Minimum value for line-search step size, must be > 0.
int InexactNewtonSolver::SetMinLinesearch(float min) {
  if (min > 0.0)  {
    LSmin = min;  
    return SUCCESS; 
  }
  else  
    return FAIL;
}

//   Use damped Newton (choice!=0) or not (choice=0)
int InexactNewtonSolver::SetDampedNewton(int choice) {
  DampNewt = choice;
  return SUCCESS; 
}

