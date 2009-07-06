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
/  Free-streaming Radiation Implicit Problem Class
/  Main Solver routine
/
/  written by: Daniel Reynolds
/  date:       June, 2009
/  modified:   
/
/  PURPOSE: Takes in relevant problem-defining parameters, as well as
/           Enzo data arrays.  These arrays may be scaled from the 
/           internal Enzo units to some that are more suitable for 
/           the implicit solve, as determined by the input parameter 
/           EScale.  Then sets up and solves the comoving, scaled 
/           equation for the free-streaming radiation, estimates the 
/           preferred time step size for the next time step, and returns 
/           the solution to the calling routine.
/
/           This routine will be called repeatedly (once per time step), 
/           so it should NOT allocate any memory.
/
************************************************************************/
#ifdef TRANSFER
#ifdef USE_HYPRE
#include "FSProb.h"
#include "CosmologyParameters.h"

/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);




int FSProb::Evolve(HierarchyEntry *ThisGrid) 
{
#ifdef USE_JBPERF
  JBPERF_START("fsprob_solve");
#endif
  
  //  if (debug)  printf("Entering FSProb::Evolve routine\n");

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;

  // declare some variables
  int i, j, k;

#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
	    MyProcessorNumber, int(MPI_id));
    return FAIL;
  }
#endif

  // start MPI timer
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  // get information from Grid
  dt = ThisGrid->GridData->ReturnTimeStep();
  told = ThisGrid->GridData->ReturnTime();
  tnew = told+dt;

  // attach radiation array to U0 vector
  float *Efold = ThisGrid->GridData->AccessRadiationFrequency0();
  U0->SetData(0, Efold);

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() == FAIL) {
    fprintf(stderr,"FSProb Solve: vector exchange_start error\n");
    return FAIL;
  }

  // get internal Enzo units (old and new time steps)
  float DenUnits, TempUnits, VelUnits, MassUnits, RadUnits;
  LenUnits0 = TimeUnits0 = LenUnits = TimeUnits = 1.0;
  DenUnits = TempUnits = VelUnits = MassUnits = RadUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits0, &TempUnits, 
	       &TimeUnits0, &VelUnits, &MassUnits, told) == FAIL) {
    fprintf(stderr,"FSProb Solve: Error in GetUnits.\n");
    return FAIL;
  }
  if (RadiationGetUnits(&RadUnits, told) == FAIL) {
    fprintf(stderr,"FSProb Solve: Error in RadiationGetUnits.\n");
    return FAIL;
  }
  EUnits0 = RadUnits*EScale;
  DenUnits = TempUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, tnew) == FAIL) {
    fprintf(stderr,"FSProb Solve: Error in GetUnits.\n");
    return FAIL;
  }
  if (RadiationGetUnits(&RadUnits, tnew) == FAIL) {
    fprintf(stderr,"FSProb Solve: Error in RadiationGetUnits.\n");
    return FAIL;
  }
  EUnits = RadUnits*EScale;

  // set a, adot, unit scalings to correct time-level values
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) == FAIL) {
      fprintf(stderr, "FSProb Solve: CosmologyComputeExpansionFactor error.\n");
      return FAIL;
    }
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) == FAIL) {
      fprintf(stderr, "FSProb Solve: CosmologyComputeExpansionFactor error.\n");
      return FAIL;
    }
    aUnits = 1.0/(1.0 + InitialRedshift);
  }

  // rescale the time to physical values
  dt *= TimeUnits;
  told *= TimeUnits0;
  tnew *= TimeUnits;
  adot /= TimeUnits;
  adot0 /= TimeUnits0;

  // begin diagnostics
  if (debug) {
    printf("\n ======================================================================\n");
    printf(" Free-Streaming Radiation Solver:\n");
  }

  // fill in the emissivity source array
  int src_set = 0;
  float *RadSrc = extsrc->GetData(0);
  if (RadSrc == NULL) {
    fprintf(stderr,"FSProb Solve: could not access Radiation source\n");
    return FAIL;
  }
#ifdef EMISSIVITY
  // if using external Emissivity field source, copy into extsrc
  if (StarMakerEmissivityField > 0) {
    // access external emissivity field 
    float *EmissivitySource = ThisGrid->GridData->AccessEmissivityField();
    if (EmissivitySource == NULL) {
      fprintf(stderr,"FSProb Solve: could not access emissivity field\n");
      return FAIL;
    }
    // copy data
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      RadSrc[i] = EmissivitySource[i];
    src_set = 1;
  }
#endif
  if (src_set == 0) {
    if (this->RadiationSource(RadSrc) != SUCCESS) {
      fprintf(stderr,"FSProb Solve: Error in RadiationSource routine\n");
      return FAIL;
    }
    src_set = 1;
  }

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() == FAIL) {
    fprintf(stderr,"FSProb Solve: vector exchange_end error\n");
    return FAIL;
  }

  // output norm of radiation sources
  float srcNorm = extsrc->rmsnorm();
  float srcMax  = extsrc->infnorm();
  if (debug) 
    printf("   emissivity norm = %g,  max = %g\n",srcNorm,srcMax);

  // rescale Enzo units with input scalings to non-dimensionalize within solver
  U0->scale_component(0,1.0/EScale);

  // output status of current solution
  float Efs_rms = U0->rmsnorm();
  float Efs_max = U0->infnorm();
  if (debug) {
    printf("   current internal (physical) values:\n");
    printf("      Efs rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   Efs_rms, Efs_rms*EUnits0, Efs_max, Efs_max*EUnits0);
  }

  // calculate initial guess at time-evolved solution
  if (this->InitialGuess(sol,U0,extsrc) != SUCCESS) {
    fprintf(stderr,"FSProb Solve: Error in InitialGuess routine\n");
    return FAIL;
  }

  // set up the linear system matrix & rhs
  float *Efnew = sol->GetData(0);
  float rhsnorm;
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, 
			Efold, RadSrc) != SUCCESS) {
    fprintf(stderr,"FSProb Solve: Error in SetupSystem routine\n");
    return FAIL;
  }
  
  // solve the free-streaming radiation problem to obtain the 
  // background propagation.
  Eint32 zed=0;
  Eint32 one=1;
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  HYPRE_StructMatrixSetBoxValues(J, ilower, iupper, stSize, entries, matentries); 

  //       assemble matrix
  HYPRE_StructMatrixAssemble(J);

  //       insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);

  //       insert sol initial guess into HYPRE vector x 
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  int xBuff, yBuff, zBuff;
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  int Zbl, Ybl, ix, iy, iz;
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = Efnew[Zbl+Ybl+xBuff+ix];
      HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
    }
  }

  //       assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

//   if (debug)  printf("Writing out matrix to file J.mat\n");
//   HYPRE_StructMatrixPrint("J.mat",J,0);

//   if (debug)  printf("Writing out rhs to file b.vec\n");
//   HYPRE_StructVectorPrint("b.vec",rhsvec,0);

//   if (debug)  printf("Writing out initial guess to file x.vec\n");
//   HYPRE_StructVectorPrint("x.vec",solvec,0);

  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_StructSolver solver;
  HYPRE_StructSolver preconditioner;
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);

  //          set preconditioner options
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
  HYPRE_StructPCGSetLogging(solver, sol_log);

  //          set solver options
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver, 
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
		      preconditioner);
  }
  else {    // ignore smg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit*500);
  }
  if (sol_tolerance != 0.0) {
    HYPRE_StructPCGSetAbsoluteTol(solver, Eflt64(sol_tolerance));
    HYPRE_StructPCGSetTol(solver, Eflt64(0.0));
  }
  HYPRE_StructPCGSetup(solver, J, rhsvec, solvec);

  //       solve the linear system
  HYPRE_StructPCGSolve(solver, J, rhsvec, solvec);

  //       extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;
  Eint32 Sits=0;
  Eint32 Pits=0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug)
    printf("   lin resid = %.1e (tol = %.1e, initial = %.1e), its = (%i,%i)\n",
	   finalresid*rhsnorm,sol_tolerance,rhsnorm,Sits,Pits);
  if (finalresid*rhsnorm > sol_tolerance) {
    fprintf(stderr," ----------------------------------------------------------------------\n");
    fprintf(stderr,"   Error: could not achieve prescribed tolerance!\n");
    fprintf(stderr," ======================================================================\n\n");
    return FAIL;
  }

  //       extract values from solution vector 
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	Efnew[Zbl+Ybl+xBuff+ix] = HYPREbuff[ix];
    }
  }

  //       destroy HYPRE solver structures
  HYPRE_StructPCGDestroy(solver);
  HYPRE_StructPFMGDestroy(preconditioner);

  // estimate the next time step size (in case it is asked)
  float diff, w, tmp;
  int x0len = LocDims[0] + GhDims[0][0] + GhDims[0][1];
  int x1len = LocDims[1] + GhDims[1][0] + GhDims[1][1];

  // enforce a solution floor on radiation values
  float Ef_floor = 1.0e-30;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
    Efnew[i] = max(Efnew[i], Ef_floor);
  
  // perform local estimates for the radiation energy relative change
  float dtfac = 0.00;   //// change this to an input parameter ////
  float dtnorm = 2.0;   //// change this to an input parameter ////
  float atol = 1.0;     //// change this to an input parameter ////
  float mindt = 0.0;    //// change this to an input parameter ////
  if (dtfac > 0.0) { 
    float loc_est = 0.0;
    if (dtnorm > 0.0) {
      for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	  for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	    w = dtfac*(sqrt(fabs(Efnew[(k*x1len + j)*x0len + i]
				*Efold[(k*x1len + j)*x0len + i])) 
			  + atol);
	    diff = Efnew[(k*x1len + j)*x0len + i] 
  	         - Efold[(k*x1len + j)*x0len + i];
	    tmp = fabs(diff/w);
	    loc_est += POW(tmp,dtnorm);
	  }
    }
    else {
      for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	  for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	    w = dtfac*(sqrt(fabs(Efnew[(k*x1len + j)*x0len + i]
				*Efold[(k*x1len + j)*x0len + i])) 
			  + atol);
	    diff = Efnew[(k*x1len + j)*x0len + i] 
  	         - Efold[(k*x1len + j)*x0len + i];
	    tmp = fabs(diff/w);
	    loc_est = (loc_est > tmp) ? loc_est : tmp;
	  }
    }
    
    // communicate to obtain overall sum/max
    float glob_est;
    int Nglobal = GlobDims[0]*GlobDims[1]*GlobDims[2];
#ifdef USE_MPI
    if (Nglobal == LocDims[0]*LocDims[1]*LocDims[2])  glob_est = loc_est;
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      if (dtnorm > 0.0) 
	MPI_Allreduce(&loc_est,&glob_est,1,DataType,MPI_SUM,MPI_COMM_WORLD);
      else
	MPI_Allreduce(&loc_est,&glob_est,1,DataType,MPI_MAX,MPI_COMM_WORLD);
    }
#else
    glob_est = loc_est;
#endif
    
    // compute overall norm
    if (dtnorm > 0.0) glob_est = POW(glob_est/Nglobal, 1.0/dtnorm);
    
    // compute time step estimate (physical units)
    dt_suggest = (glob_est == 0.0) ? huge_number : dt/glob_est;
    dt_suggest = min(dt_suggest, huge_number);
    
    // limit maximum growth per step
    dt_suggest = min(dt_suggest, 1.1*dt);
    
    // rescale dt estimates to normalized values
    dt_suggest /= TimeUnits;
    
    // account for min/max time step size (according to user)
    dt_suggest = max(dt_suggest, mindt);
    dt_suggest = min(dt_suggest, maxdt);
    
    if (debug)  printf("   dt_suggest = %8.2e\n",dt_suggest);
    ThisGrid->GridData->SetMaxRadiationDt(dt_suggest);

  } // end if (dtfac > 0.0)
  else if (maxdt > 0.0) {
    dt_suggest = maxdt;
    if (debug)  printf("   dt_suggest = %8.2e\n",dt_suggest);
    ThisGrid->GridData->SetMaxRadiationDt(dt_suggest);
    
  }

  // rescale results back to Enzo units
  sol->scale_component(0,EScale);
  
  // rescale the time to normalized values
  dt /= TimeUnits;
  told /= TimeUnits0;
  tnew /= TimeUnits;
  adot *= TimeUnits;
  adot0 *= TimeUnits0;

  // update Enzo data with new values (Enzo pointer to Ef is in U0)
  U0->copy(sol);

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  FStime += ftime-stime;
  if (debug)  printf("   FS cumulative wall time = %g\n",FStime);

  if (debug)
    printf(" ======================================================================\n\n");

#ifdef USE_JBPERF
  JBPERF_STOP("fsprob_solve");
#endif

  // return success
  return SUCCESS;
}
#endif
#endif
