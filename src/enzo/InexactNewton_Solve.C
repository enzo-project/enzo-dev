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
/  Inexact Newton Solve Routine
/
/  written by: Daniel R. Reynolds
/  date:       March, 2006
/  modified1:
/
/  PURPOSE: Implements a damped inexact Newton method:  Given an initial 
/           guess x0, while nonlinear residual is too large, do
/                  (i) solve J(x)*s = -f(x) + r
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
/                  6  ->  ||f(x)||_rms / ||x0||_rms < Ntol
/                  7  ->  ||f(x)||_inf / ||x0||_inf < Ntol
/           where ||.||_rms denotes the root-mean-squared norm (2-norm 
/           divided by number of values) and ||.||_inf denotes the 
/           infinity norm (or maximum absolute value).  Choices 1 and 2 
/           require raw residual values below Ntol; choices 2 and 3 
/           require relative residual decrease below Ntol; choices 4 
/           and 5 weight the residual by the magnitude of the solution; 
/           choices 6 and 7 weight the residual by the magnitude of the 
/           previous time-level solution.  Control over the test type 
/           and Ntol are given through the SetNewtonNorm() and 
/           SetNewtonTolerance() routines.
/
/           The inexactness parameter, r, in step (i) is given through 
/           the InexactNewtonForce() routine.  The linesearch step size
/           l is found in the LinesearchStepSize() routine.  The linear 
/           solver used in step (i) is provided by the problem-specific 
/           lsolve() routine.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "performance.h"
#include "InexactNewton.h"
#include "global_data.h"


int InexactNewtonSolver::Solve(NonlinearProblemABC *prob, EnzoVector *x)
{
//   if (debug)  printf("Entering InexactNewtonSolver::Solve routine\n");

  // local variable (for convergence test)
  float normscale, fnormtest;

  // get initial nonlinear residual and norm
  if (prob->nlresid(fvec, x) != SUCCESS) {
    fprintf(stderr,"Error in problem nlresid routine.\n");
     return FAIL;
  }
  fnorm = fnorm0 = fvec->rmsnorm();

  // set convergence norm scalings if needed
  if      (NtolNorm == 2)  normscale = fnorm0;
  else if (NtolNorm == 3)  normscale = fvec->infnorm();
  else if (NtolNorm == 6)  normscale = x->rmsnorm();
  else if (NtolNorm == 7)  normscale = x->infnorm();
  if (normscale < 1e-16)   normscale = 1.0;   // just in case


  // base convergence test scaling off of NtolNorm
  if      (NtolNorm == 0)  fnormtest = fnorm;
  else if (NtolNorm == 1)  fnormtest = fvec->infnorm();
  else if (NtolNorm == 2)  fnormtest = fnorm/normscale;
  else if (NtolNorm == 3)  fnormtest = fvec->infnorm()/normscale;
  else if (NtolNorm == 4)  fnormtest = fnorm/(x->rmsnorm());
  else if (NtolNorm == 5)  fnormtest = (fvec->infnorm())/(x->infnorm());
  else if (NtolNorm == 6)  fnormtest = fnorm/normscale;
  else if (NtolNorm == 7)  fnormtest = (fvec->infnorm())/normscale;

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // begin Newton iterations (if needed)
  Niters = 0;
  if (fnormtest >= Ntol) {
    for (Niters=0; Niters<Nmax_it; Niters++) {
    
      if (debug) {
	printf("  -----------------------------------------------------------------------\n");	
	printf("   Newton iteration %"ISYM", ||f|| = %g (tol = %g)\n",
	       Niters,fnormtest,Ntol);
	printf("  -----------------------------------------------------------------------\n");	
      }
      
      // obtain inexact Newton force
      if (InexactNewtonForce() != SUCCESS) {
	fprintf(stderr, "Error: Inexact Newton Force error\n");
	return FAIL;
      }
      
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      // set up Jacobian for Newton solve
      if (prob->lsetup(x) != SUCCESS) {
	fprintf(stderr,"Error: Newton linear system setup error\n");
	return FAIL;
      }
      
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      // solve problem-specific inexact Newton system, J(u)*dx=f
      if (prob->lsolve(dx, fvec, x, INforce) != SUCCESS) {
	fprintf(stderr,"Error: Newton linear solver error\n");
	return FAIL;
      }
      
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      if (DampNewt != 0) {
	// perform linesearch globalization
	if (LinesearchStepSize(prob, x, dx) != SUCCESS) {
// 	  fprintf(stderr,"Error: Line-search step size error\n");
// 	  return FAIL;
	}
	// NOTE: after a successful linesearch globalization, the temporary 
	// vector tmp1 contains the Newton updated state x = x - LSsize*dx; 
	// fvec contains the corresponding nonlinear residual, and fnorm 
	// contains the RMS norm of fvec.

	// copy vtmp1 to x
	if (x->copy(tmp1) != SUCCESS) {
	  fprintf(stderr,"Error in EnzoVector copy routine\n");
	  return FAIL;
	}
      }
      else {
	// perform the update manually
	LSsize = 1.0;
	if (x->axpy(-LSsize, dx) != SUCCESS) {
	  fprintf(stderr, "Error in EnzoVector axpy routine!\n");
	  return FAIL;
	}
	if (prob->nlresid(fvec, x) != SUCCESS) {
	  fprintf(stderr, "Error in Problem nlresid routine!\n");
	  return FAIL;
	}
	fnorm = fvec->rmsnorm();
      }

      // base convergence test scaling off of NtolNorm
      if      (NtolNorm == 0)  fnormtest = fnorm;
      else if (NtolNorm == 1)  fnormtest = fvec->infnorm();
      else if (NtolNorm == 2)  fnormtest = fnorm/normscale;
      else if (NtolNorm == 3)  fnormtest = fvec->infnorm()/normscale;
      else if (NtolNorm == 4)  fnormtest = fnorm/(x->rmsnorm());
      else if (NtolNorm == 5)  fnormtest = (fvec->infnorm())/(x->infnorm());
      else if (NtolNorm == 6)  fnormtest = fnorm/normscale;
      else if (NtolNorm == 7)  fnormtest = (fvec->infnorm())/normscale;
      
      // check for convergence, otherwise set old residual value for next pass
      if (fnormtest < Ntol) {
	Niters += 1;
	break;
      }
      else fnorm0 = fnorm;
    }
  }
  
  // output final residual and return
  if (debug) {
    printf("  -----------------------------------------------------------------------\n");
    printf("   Final Newton residual: %i iters, ||f|| = %g (tol = %g)\n", 
	   Eint32(Niters), fnormtest, Ntol);
    printf("  =======================================================================\n");
  }

  return SUCCESS;
}
