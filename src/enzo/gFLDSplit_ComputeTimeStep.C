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
/  Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Computes the rad-hydro time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
/
/           Also, the argument 'flag' determines which entries we use to 
/           compute the time step:
/                0 -> use radiation only
/                1 -> use chemistry + energy correction only
/                2 -> use both
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

 
 

float gFLDSplit::ComputeTimeStep(EnzoVector *uold, EnzoVector *unew, int flag)
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

  // set internal dtfactor array based on user inputs and flag:
  float dtfactor[Nchem+2];
  int i;
  for (i=0; i<2+Nchem; i++)  dtfactor[i]=huge_number;
  if (flag == 0)  
    dtfactor[0] = dtfac[0];
  if (flag == 1)  {
    dtfactor[1]=dtfac[1];
    for (i=2; i<2+Nchem; i++)  dtfactor[i]=dtfac[2];
  }
  if (flag == 2) {
    for (i=0; i<2; i++)        dtfactor[i]=dtfac[i];
    for (i=2; i<2+Nchem; i++)  dtfactor[i]=dtfac[2];
  }


  // Set time step depending on how it has been set up by the user:
  //    If dtfactor is set for any species, compute maximum time step 
  //    as estimate allowing dtfactor relative change.  This relative 
  //    change is estimated as follows:
  //       dt_new = dt_old / relerr_fac
  //    where relerr_fac gives the ratio between an estimated 
  //    local truncation error and the desired relative change:
  //       relerr_fac = || (unew - uold) / w ||_p
  //    with the scaling vector w given by
  //       w = dtfactor*[sqrt(|unew*uold|) + atol]
  //    and where we have the following parameters:
  //       p - norm choice (input), 0->max norm, otherwise the p-norm
  //           **all p-norms here divide by the total number of cells**
  //       dtfactor - desired relative change per step (input)
  //       atol - 1e-3 (assumes units are all normalized)
  //    For the gas energy correction, this is different since we do 
  //    not have uold: 
  //       relerr_fac = || unew / w ||_p
  //       w = dtfactor*(|unew| + atol).
  //
  //    If dtfactor is not set for any species, use a heuristic to shoot 
  //    for 2 to 3 newton iterations per time step, with linesearch 
  //    step length equal to 1.
  float dt_est = huge_number;    // max time step (normalized units)
  float test = dtfactor[0];
  for (i=0; i<2+Nchem; i++)  test = min(dtfactor[i],test);
  if (test != huge_number) {

    // initialize variables
    float diff, w, tmp, atol;
    int j, k, l;
    int x0len = Nx + ghXl + ghXr;
    int x1len = Ny + ghYl + ghYr;
    float loc_est[Nvar];

    // perform local estimates for the radiation energy relative change
    loc_est[0] = 0.0;
    if (dtfactor[0] != huge_number) {
      float *Eold = uold->GetData(0);
      float *Enew = unew->GetData(0);
      atol = 0.001; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfactor[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
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
	      w = dtfactor[0]*(sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
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
    if (dtfactor[1] != huge_number) {
      float *ec = unew->GetData(1);
      atol = 0.001; // assumes values are normalized
      if (dtnorm > 0.0) {
	for (k=ghZl; k<Nz+ghZl; k++) 
	  for (j=ghYl; j<Ny+ghYl; j++)
	    for (i=ghXl; i<Nx+ghXl; i++) {
	      w = dtfactor[1]*(fabs(ec[(k*x1len + j)*x0len + i]
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
	      w = dtfactor[1]*(fabs(ec[(k*x1len + j)*x0len + i]
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
      if (dtfactor[l] != huge_number) {
	niold = uold->GetData(l);
	ninew = unew->GetData(l);
	atol = 0.001; // assumes values are normalized
	if (dtnorm > 0.0) {
	  for (k=ghZl; k<Nz+ghZl; k++) 
	    for (j=ghYl; j<Ny+ghYl; j++)
	      for (i=ghXl; i<Nx+ghXl; i++) {
		w = dtfactor[l]*(sqrt(fabs(ninew[(k*x1len + j)*x0len + i]
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
		w = dtfactor[l]*(sqrt(fabs(ninew[(k*x1len + j)*x0len + i]
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
      printf("  gFLDSplit_ComputeTimestep: (E, e, ni) dt_est = (");
      for (l=0; l<Nvar; l++) {
	if (dt_est_var[l] == huge_number/TimeUnits)
	  printf(" -------- ");
	else  
	  printf(" %8.2e ",dt_est_var[l]);
      }
      printf(")\n");
    }
  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;

}

#endif
