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
/  Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       March, 2007
/  modified1:  
/
/  PURPOSE: Computes the rad-hydro time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called just before gFLDProblem::Return, so 
/           it sees the same units as the rest of the solver module.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"

 
 

float gFLDProblem::ComputeTimeStep(EnzoVector *uold, EnzoVector *unew, 
				   int NewtIts, float FStep, float FResid)
{
  // get local mesh description
  int Nx, Ny, Nz, Nvar, ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  unew->size(&Nx, &Ny, &Nz, &Nvar, &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (Nx != LocDims[0]) 
    ENZO_FAIL("ComputeTimeStep error: x0 vector dims do not match");
  if (Ny != LocDims[1]) 
    ENZO_FAIL("ComputeTimeStep error: x1 vector dims do not match");
  if (Nz != LocDims[2]) 
    ENZO_FAIL("ComputeTimeStep error: x2 vector dims do not match");
  if (Nvar != (2+Nchem)) 
    ENZO_FAIL("ComputeTimeStep error: nspecies dims do not match");
  if ((Nx+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("ComputeTimeStep error: x0 vector sizes do not match");
  if ((Ny+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("ComputeTimeStep error: x1 vector sizes do not match");
  if ((Nz+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("ComputeTimeStep error: x2 vector sizes do not match");


  // Set time step depending on how it has been set up by the user:
  //    If dtfac is set for any species, compute maximum time step 
  //    as estimate allowing dtfac relative change.  This relative 
  //    change is estimated as follows:
  //       dt_new = dt_old / relerr_fac
  //    where relerr_fac gives the ratio between an estimated 
  //    local truncation error and the desired relative change:
  //       relerr_fac = || (unew - uold) / w ||_p
  //    with the scaling vector w given by
  //       w = dtfac*[sqrt(|unew*uold|) + atol]
  //    and where we have the following parameters:
  //       p - norm choice (input), 0->max norm, otherwise the p-norm
  //           **all p-norms here divide by the total number of cells**
  //       dtfac - desired relative change per step (input)
  //       atol - 1e-3 (assumes units are all normalized)
  //    For the gas energy correction, this is different since we do 
  //    not have uold: 
  //       relerr_fac = || unew / w ||_p
  //       w = dtfac*(|unew| + atol).
  //
  //    If dtfac is not set for any species, use a heuristic to shoot 
  //    for 2 to 3 newton iterations per time step, with linesearch 
  //    step length equal to 1.
  float dt_est = huge_number;
  float test = dtfac[0];
  int i;
  for (i=0; i<2+Nchem; i++)  test = min(dtfac[i],test);
  if (test != huge_number) {

    // initialize variables
    float diff, w, tmp, atol;
    int i, j, k, l;
    int x0len = Nx + ghXl + ghXr;
    int x1len = Ny + ghYl + ghYr;
    float loc_est[Nvar];

    // perform local estimates for the radiation energy relative change
    loc_est[0] = 0.0;
    if (dtfac[0] != huge_number) {
      float *Eold = uold->GetData(0);
      float *Enew = unew->GetData(0);
      atol = 0.001; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[0] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
				     *Eold[(k*x1len + j)*x0len + i])) 
			    + atol);
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[0] = (loc_est[0] > tmp) ? loc_est[0] : tmp;
	    }
      }
    }

    // perform estimates for the gas energy
    loc_est[1] = 0.0;
    if (dtfac[1] != huge_number) {
      float *ec = unew->GetData(1);
      atol = 0.001; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[1]*(fabs(ec[(k*x1len + j)*x0len + i]
			       + eh[(k*x1len + j)*x0len + i]/ecScale) 
			    + atol);
	      diff = ec[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[1] += POW(tmp,dtnorm);
	    }
      }
      else {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfac[1]*(fabs(ec[(k*x1len + j)*x0len + i]
			       + eh[(k*x1len + j)*x0len + i]/ecScale) 
			    + atol);
	      diff = ec[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      loc_est[1] = (loc_est[1] > tmp) ? loc_est[1] : tmp;
	    }
      }
    }

    // perform estimates for the chemistry
    float *niold, *ninew;
    for (l=2; l<=Nchem+1; l++) {
      loc_est[l] = 0.0;
      if (dtfac[2] != huge_number) {
	niold = uold->GetData(l);
	ninew = unew->GetData(l);
	atol = 0.001; // assumes values are normalized
	if (dtnorm > 0.0) {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		w = dtfac[2]*(sqrt(fabs(ninew[(k*x1len + j)*x0len + i]
				       *niold[(k*x1len + j)*x0len + i])) 
				+ atol);
		diff = ninew[(k*x1len + j)*x0len + i] 
		     - niold[(k*x1len + j)*x0len + i];
		tmp = fabs(diff/w);
		loc_est[l] += POW(tmp,dtnorm);
	      }
	}
	else {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		w = dtfac[2]*(sqrt(fabs(ninew[(k*x1len + j)*x0len + i]
				       *niold[(k*x1len + j)*x0len + i])) 
				+ atol);
		diff = ninew[(k*x1len + j)*x0len + i]
   	 	     - niold[(k*x1len + j)*x0len + i];
		tmp = fabs(diff/w);
		loc_est[l] = (loc_est[l] > tmp) ? loc_est[l] : tmp;
	      }
	}
      }
    }

    // communicate to obtain overall sum/max
    float glob_est[Nvar];
    int Nglobal = GlobDims[0]*GlobDims[1]*GlobDims[2];
#ifdef USE_MPI
    if (Nglobal == Nx*Ny*Nz) 
      for (l=0; l<Nvar; l++)  glob_est[l] = loc_est[l];
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg vars = Nvar;
      if (dtnorm > 0.0) 
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_SUM,MPI_COMM_WORLD);
      else
	MPI_Allreduce(&loc_est,&glob_est,vars,DataType,MPI_MAX,MPI_COMM_WORLD);
    }
#else
    for (l=0; l<Nvar; l++)  glob_est[l] = loc_est[l];
#endif

    // compute overall norms
    if (dtnorm > 0.0) 
      for (l=0; l<Nvar; l++) {
	glob_est[l] /= Nglobal;
	glob_est[l] = POW(glob_est[l],1.0/dtnorm);
      }

    // compute variable-specific time step estimates (physical units)
    float dt_est_var[Nvar];
    for (l=0; l<Nvar; l++) {
      dt_est_var[l] = (glob_est[l] == 0.0) ? huge_number : dt/glob_est[l];
      dt_est_var[l] = min(dt_est_var[l], huge_number);
    }

    // set estimated time step as minimum of component time steps
    dt_est = maxdt*TimeUnits;    // max time step estimate (physical units)
    for (l=0; l<Nvar; l++) {
      dt_est = min(dt_est, dt_est_var[l]);
    }

    // limit maximum growth per step
    dt_est = min(dt_est, 1.1*dt);    // time step growth (physical units)

    // rescale dt estimates to normalized values
    dt_est /= TimeUnits;
    for (l=0; l<Nvar; l++)  dt_est_var[l] /= TimeUnits;

    // account for min/max time step size (according to user)
    dt_est = max(dt_est, mindt);
    dt_est = min(dt_est, maxdt);

    if (debug) {
      printf("gFLDProblem_ComputeTimestep: (E, e, ni) dt_est = (");
      for (l=0; l<Nvar; l++)  printf(" %8.2e ",dt_est_var[l]);
      printf(")\ngFLDProblem_ComputeTimeStep: dt_est = %g\n",dt_est);
    }
  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;
}

#endif
