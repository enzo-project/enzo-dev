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
/  Evolve Routine for split radiation-chemistry system
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "gFLDSplit.h"
#include "InexactNewton.h"
#include "CosmologyParameters.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);




int gFLDSplit::Evolve(HierarchyEntry *ThisGrid, float deltat)
{

//   if (debug)  printf("Entering gFLDSplit::Evolve routine\n");

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;

#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    fprintf(stderr, "ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", 
	    MyProcessorNumber, int(MPI_id));
    ENZO_FAIL("Error in gFLDSplit_Evolve");
  }
#endif

  // in case MPI is not included
#ifndef MPI_INT
  int MPI_COMM_WORLD = 0;
#endif

  // start MPI timer
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  // Set pointers to each variable (zero out fluid energy correction 
  // since that is not internal to Enzo)
  float *RadiationEnergy = NULL;
  int i;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    FluidEnergyCorrection[i] = 0.0;
  vx = ThisGrid->GridData->AccessVelocity1();
  if (vx == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not obtain velocity1");
  vy = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not obtain velocity2");
  vz = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not obtain velocity3");
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not obtain density");
  if (DualEnergyFormalism) {
    eh = ThisGrid->GridData->AccessGasEnergy();
    if (eh == NULL) 
      ENZO_FAIL("gFLDSplit Evolve: could not obtain fluid energy");
  }
  else {
    eh = ThisGrid->GridData->AccessTotalEnergy();
    if (eh == NULL) 
      ENZO_FAIL("gFLDSplit Evolve: could not obtain fluid energy");
  }
  RadiationEnergy = ThisGrid->GridData->AccessRadiationFrequency0();
  if (RadiationEnergy == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not obtain Radiation energy");
  // "access" all chemical species (some will be NULL); this helps for 
  // problems in which we only do Hydrogen chemistry internally to this 
  // module, but Helium species are present.
  float *nHI    = NULL;
  float *nHII   = NULL;
  float *nHeI   = NULL;
  float *nHeII  = NULL;
  float *nHeIII = NULL;
  float *ne     = NULL;
  nHI    = ThisGrid->GridData->AccessHIDensity();
  nHII   = ThisGrid->GridData->AccessHIIDensity();
  nHeI   = ThisGrid->GridData->AccessHeIDensity();
  nHeII  = ThisGrid->GridData->AccessHeIIDensity();
  nHeIII = ThisGrid->GridData->AccessHeIIIDensity();
  ne     = ThisGrid->GridData->AccessElectronDensity();
  //    check that we accessed the required species for Nchem
  if (Nchem > 0) 
    if (nHI == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HI density!");
  if (Nchem > 1) {
    if (nHeI == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HeI density!");
    if (nHeII == NULL) 
      ENZO_FAIL("EvolveRadHydro error: cannot access HeII density!");
  }

  // Get time-related information
  dt = deltat;
  told = ThisGrid->GridData->ReturnTime();
  tnew = told+dt;

  // get internal Enzo units (old time step)
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0=LenUnits0=TempUnits=TimeUnits=VelUnits=MassUnits=1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, told) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) == FAIL) 
    ENZO_FAIL("Error in RadiationGetUnits.");
  // incorporate Enzo units with implicit solver unit scaling
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits0 = RadUnits*ErScale;
  NiUnits0 = (Nchem == 0) ? NiScale : DenUnits0/mp*NiScale;

  // set a, adot, unit scalings to correct time-level values
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");

  // get/store internal Enzo units (new time step)
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, tnew) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) == FAIL) 
    ENZO_FAIL("Error in RadiationGetUnits.");
  // incorporate Enzo units with implicit solver unit scaling
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits*ecScale;
  NiUnits = (Nchem == 0) ? NiScale : DenUnits/mp*NiScale;

  // set a, adot to correct time-level values
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");

  // initialize external sources to 0
  extsrc->constant(0.0);

  // access/fill radiation source array 
  float *RadSrc = extsrc->GetData(0);
  if (RadSrc == NULL) 
    ENZO_FAIL("gFLDSplit Evolve: could not access Radiaton source");
  int eta_set = 0;
#ifdef EMISSIVITY
  // if using external Emissivity field source, copy into extsrc
  if (StarMakerEmissivityField > 0) {
    // access external emissivity field 
    float *EmissivitySource = ThisGrid->GridData->AccessEmissivity0();
    if (EmissivitySource == NULL) 
      ENZO_FAIL("gFLDSplit Evolve: could not access emissivity field");
    // copy data
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      RadSrc[i] = EmissivitySource[i];

    eta_set = 1;
  }
#endif
  //   compute emissivity (if not yet set)
  if (eta_set == 0) {
    if (this->RadiationSource(RadSrc, &tnew) != SUCCESS) 
      ENZO_FAIL("LocRHS: Error in RadiationSource routine");
  }
  float srcNorm = extsrc->rmsnorm_component(0);
  float srcMax  = extsrc->infnorm_component(0);
  if (debug) {
    printf("\n gFLDSplit Evolve:\n");
    printf("   emissivity norm = %g,  max = %g\n",srcNorm,srcMax);
  }


  // adjust chemical species HI, HeI, HeII to correspond to rho (since 
  // rho and ni are advected differently in hydro solver for some reason)
  float rhochem;
  // Hydrogen chemistry 
  if (Nchem > 0) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      // first ensure that no densities are negative
      nHI[i]  = max(0.0, nHI[i]);
      nHII[i] = max(0.0, nHII[i]);

      // set rhochem as the total H 'density' according to H* species
      rhochem = nHI[i] + nHII[i];

      // update HI as appropriate fraction of 'true' density
      nHI[i] *= rho[i]*HFrac/rhochem;

      // correct if precision isn't good enough, and all rho is HI
      nHI[i]  = min(nHI[i],rho[i]*HFrac);
      nHII[i] = max(0.0, rho[i]*HFrac - nHI[i]);
    }
  }

  // Helium chemistry 
  // (or if we are running Hydrogen chemistry, but Helium is present in code)
  if ((Nchem > 1) || (HFrac < 1.0)) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {

      // first ensure that no densities are negative
      nHeI[i]   = max(0.0, nHeI[i]);
      nHeII[i]  = max(0.0, nHeII[i]);
      nHeIII[i] = max(0.0, nHeIII[i]);

      // set rhochem as the total He 'density' according to He* species
      rhochem = nHeI[i] + nHeII[i] + nHeIII[i];

      // update HeI, HeII as appropriate fractions of 'true' density
      nHeI[i]   *= rho[i]*(1.0-HFrac)/rhochem;
      nHeII[i]  *= rho[i]*(1.0-HFrac)/rhochem;
      nHeIII[i] *= max(0.0, rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i]);
    }
  }

  // attach arrays to U0 vector
  U0->SetData(0, RadiationEnergy);
  U0->SetData(1, FluidEnergyCorrection);
  if (Nchem > 0)  U0->SetData(2, nHI);
  if (Nchem > 1) {
    U0->SetData(3, nHeI);
    U0->SetData(4, nHeII);
  }

  // rescale Enzo units with input scalings to non-dimensionalize within solver
  U0->scale_component(0,1.0/ErScale);
  U0->scale_component(1,1.0/ecScale);
  int ns;
  for (ns=1; ns<=Nchem; ns++)  U0->scale_component(ns+1, 1.0/NiScale);

  // have U0 begin communication of neighbor information
  if (U0->exchange_start() == FAIL) 
    ENZO_FAIL("gFLDSplit Evolve: vector exchange_start error");

  // output typical/maximum values
  float UTypVals[Nchem+2];
  float UMaxVals[Nchem+2];
  UTypVals[0] = U0->rmsnorm_component(0);
  UMaxVals[0] = U0->infnorm_component(0);
  for (ns=1; ns<=Nchem; ns++) {
    UTypVals[1+ns] = U0->rmsnorm_component(ns+1);
    UMaxVals[1+ns] = U0->infnorm_component(ns+1);
  }

  //    set fluid energy correction "typical" value (since ec0=0)
  UTypVals[1] = 0.0;  UMaxVals[1] = 0.0;
  float dtmp;
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
    UTypVals[1] += eh[i]*eh[i];
    UMaxVals[1] = max(UMaxVals[1],eh[i]*eh[i]);
  }
  UTypVals[1] = sqrt(UTypVals[1]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
  UMaxVals[1] = sqrt(UMaxVals[1]);
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg one = 1;
  MPI_Allreduce(&(UMaxVals[1]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  UMaxVals[1] = dtmp;
  MPI_Allreduce(&(UTypVals[1]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  UTypVals[1] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
#endif
  UTypVals[1] /= ecScale;
  UMaxVals[1] /= ecScale;

  if (debug) {
    printf("   current internal (physical) quantities:\n");
    printf("      Eg rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[0],UTypVals[0]*ErUnits0, UMaxVals[0],UMaxVals[0]*ErUnits0);
    printf("      ec rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[1],UTypVals[1]*ecUnits, UMaxVals[1],UMaxVals[1]*ecUnits);
    if (Nchem == 1) {
      printf("     nHI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[2],UTypVals[2]*NiUnits0, UMaxVals[2],UMaxVals[2]*NiUnits0);
    }
    if (Nchem == 3) {
      printf("    nHeI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[3],UTypVals[3]*NiUnits0, UMaxVals[3],UMaxVals[3]*NiUnits0);
      printf("   nHeII rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[4],UTypVals[4]*NiUnits0, UMaxVals[4],UMaxVals[4]*NiUnits0);
    }
  }


//   int buff = (GhDims[2][0]*ArrDims[1] + GhDims[1][0])*ArrDims[0] + GhDims[0][0];
//   if (debug) {
//     printf("  columns of values (1st 4 of each):\n");
//     printf("      Eg:  %8.2e  %8.2e  %8.2e  %8.2e\n",
// 	   RadiationEnergy[buff], RadiationEnergy[buff+1], 
// 	   RadiationEnergy[buff+2], RadiationEnergy[buff+3]);
//     printf("      eh:  %8.2e  %8.2e  %8.2e  %8.2e\n",
// 	   eh[buff], eh[buff+1], eh[buff+2], eh[buff+3]);
//     printf("     nHI:  %8.2e  %8.2e  %8.2e  %8.2e\n",
// 	   nHI[buff], nHI[buff+1], nHI[buff+2], nHI[buff+3]);
//   }


  // rescale dt, told, tnew, adot to physical values
  dt    *= TimeUnits;
  told  *= TimeUnits;
  tnew  *= TimeUnits;
  adot  /= TimeUnits;
  adot0 /= TimeUnits;

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() == FAIL) 
    ENZO_FAIL("gFLDSplit Evolve: vector exchange_end error");



  ////////////////////////////////////
  // Problem Solve Phase

  //   enforce boundary conditions on old time step vector
  if (this->EnforceBoundary(U0) == FAIL) 
    ENZO_FAIL("ERROR: EnforceBoundary failure!!");

  //   obtain initial guess for time-evolved solution
  if (this->InitialGuess(sol) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("ERROR: InitialGuess failure!!");
  }

  // do not use initial guess for radiation equation
  sol->copy_component(U0,0);

  //   enforce boundary conditions on new time step initial guess vector
  if (this->EnforceBoundary(sol) == FAIL) 
    ENZO_FAIL("ERROR: EnforceBoundary failure!!");

//   // write out radiation field initial guess
//   if (debug)  printf("Writing out initial guess to file guess.vec\n");
//   sol->writeall("guess.vec",0);

//   // write out gas energy initial guess
//   if (debug)  printf("Writing out energy initial guess to file energy.vec\n");
//   sol->writeall("energy.vec",1);


  //   compute updated opacities
  if (this->Opacity(OpacityE, &tnew, sol) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in Opacity routine");

//   if (debug)
//     printf("      kE:  %8.2e  %8.2e  %8.2e  %8.2e\n",
// 	   OpacityE[buff]*NiUnits, OpacityE[buff+1]*NiUnits, 
// 	   OpacityE[buff+2]*NiUnits, OpacityE[buff+3]*NiUnits);

  //   compute the gas 'Temperature' at old and new time steps
  if (this->ComputeTemperature(Temperature0, U0) == FAIL) 
    ENZO_FAIL("LocRHS: Error in ComputeTemperature routine");
  if (this->ComputeTemperature(Temperature, sol) == FAIL) 
    ENZO_FAIL("LocRHS: Error in ComputeTemperature routine");

//   // output typical/maximum temperature values
//   UTypVals[1] = 0.0;  UMaxVals[1] = 0.0;
//   for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
//     UTypVals[1] += Temperature[i]*Temperature[i];
//     UMaxVals[1] = max(UMaxVals[1],Temperature[i]*Temperature[i]);
//   }
//   UTypVals[1] = sqrt(UTypVals[1]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
//   UMaxVals[1] = sqrt(UMaxVals[1]);
// #ifdef USE_MPI
//   MPI_Allreduce(&(UMaxVals[1]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
//   UMaxVals[1] = dtmp;
//   MPI_Allreduce(&(UTypVals[1]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
//   UTypVals[1] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
// #endif
//   if (debug) 
//     printf("      Temp rms = %13.7e, max = %13.7e\n",UTypVals[1], UMaxVals[1]);
  
  

  //    access updated radiation energy array
  float *Eg_new = sol->GetData(0);


#ifdef USE_HYPRE

  // set up and solve radiation equation
  float rhsnorm;
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, RadiationEnergy, Eg_new, 
			OpacityE, Temperature, Temperature0, RadSrc) != SUCCESS) 
    ENZO_FAIL("FSProb Solve: Error in SetupSystem routine");
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 
    
  //       assemble matrix
  HYPRE_StructMatrixAssemble(P);
    
  //       insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);

  //       set linear solver tolerance (rescale to relative residual and not actual)
  Eflt64 delta = sol_tolerance / rhsnorm;
  delta = min(delta, 1.0e-6);

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
	HYPREbuff[ix] = 0.0;
      HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
    }
  }
    
  //       assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);
    
//   if (debug)  printf("Writing out matrix to file P.mat\n");
//   HYPRE_StructMatrixPrint("P.mat",P,1);

//   if (debug)  printf("Writing out rhs to file b.vec\n");
//   HYPRE_StructVectorPrint("b.vec",rhsvec,1);

//   if (debug)  printf("Writing out initial guess to file x.vec\n");
//   HYPRE_StructVectorPrint("x.vec",solvec,1);

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
    
  //          set solver options
  HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
  HYPRE_StructPCGSetLogging(solver, sol_log);
  HYPRE_StructPCGSetRelChange(solver, 1);
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
  if (delta != 0.0)   HYPRE_StructPCGSetTol(solver, delta);
  HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
    
  //       solve the linear system
  if (debug)
    printf(" ----------------------------------------------------------------------\n");
  HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);
    
//   if (debug)  printf("Writing out solution to file s.vec\n");
//   HYPRE_StructVectorPrint("s.vec",solvec,0);

  //       extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;
  Eint32 Sits=0;
  Eint32 Pits=0;
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver, &Sits);
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug) 
    printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
	   finalresid*rhsnorm, sol_tolerance, rhsnorm, Sits, Pits);
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
    // if the final residual is too large, or is nan, quit
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {
      fprintf(stderr,"   Error: could not achieve prescribed tolerance!\n");

      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P.mat\n");
      HYPRE_StructMatrixPrint("P.mat",P,0);

      if (debug)  printf("Writing out rhs to file b.vec\n");
      HYPRE_StructVectorPrint("b.vec",rhsvec,0);

      if (debug)  printf("Writing out current solution to file x.vec\n");
      HYPRE_StructVectorPrint("x.vec",solvec,0);

      // dump module parameters to disk
      this->Dump(sol);

      fprintf(stderr," ======================================================================\n\n");
      ENZO_FAIL("Error in gFLDSplit_Evolve");
    }
  }
  if (debug)	printf(" ======================================================================\n\n");

    
  //       extract values from solution vector, adding them to solution
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
 	Eg_new[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
    }
  }
    
  //       destroy HYPRE solver structures
  HYPRE_StructPCGDestroy(solver);
  HYPRE_StructPFMGDestroy(preconditioner);
  
#else  // ifdef USE_HYPRE

  ENZO_FAIL("gFLDSplit_Evolve ERROR: module requires USE_HYPRE to be set!");
  
#endif


//   // write out updated radiation field to disk
//   if (debug)  printf("Writing out solution to file sol.vec\n");
//   sol->writeall("sol.vec",0);


  // enforce a solution floor on radiation
  //   first determine minimum floating-point roundoff
  float epsilon=1.0;
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;   // smallest representable number
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    Eg_new[i] = max(Eg_new[i],epsilon);


  // subcycle the chemistry and gas energy equations
  float thisdt, dtchem2;
  float tchem = told;
  float *ecsrc   = extsrc->GetData(1);
  float *HIsrc   = NULL;
  float *HeIsrc  = NULL;
  float *HeIIsrc = NULL;
  if (Nchem > 0)  HIsrc = extsrc->GetData(2);
  if (Nchem > 1) {
    HeIsrc = extsrc->GetData(3);
    HeIIsrc = extsrc->GetData(4);
  }
  float *Opacity_new = Temperature;
  float *Opacity_old = OpacityE;
  float *sol_ec   = sol->GetData(1);
  float *sol_HI   = NULL;
  float *sol_HeI  = NULL;
  float *sol_HeII = NULL;
  if (Nchem > 0)  sol_HI = sol->GetData(2);
  if (Nchem > 1) {
    sol_HeI  = sol->GetData(3);
    sol_HeII = sol->GetData(4);
  }
  float *tmp, factor, epsilon2, *eh_tot, *eh_gas;
  for (int chemstep=0; chemstep<=100; chemstep++) {

    // update tchem
    thisdt = min(dtchem, dt);          // do not exceed radiation dt
    thisdt = max(thisdt, dt/100);      // take at most 100 steps
    tchem += thisdt;                   // update chemistry time
    //    check that we don't exceed radiation time
    if (tchem >= tnew) {
      thisdt = tnew - (tchem - thisdt);  // get max time step
      tchem = tnew;                      // set updated time
    }
    if (debug) {
      printf("  subcycled chem %"ISYM": dt=%7.1e, t=%7.1e (rad dt=%7.1e, t=%7.1e)\n",chemstep,thisdt/TimeUnits,tchem/TimeUnits,dt/TimeUnits,tnew/TimeUnits);
    }


    //   fill in the gas energy and chemistry source terms
    if (this->GasEnergySource(ecsrc, &tchem) != SUCCESS)
      ENZO_FAIL("gFLDSplit Evolve: Error in GasEnergySource routine");
    if (this->ChemistrySource(HIsrc, HeIsrc, HeIIsrc, &tchem) != SUCCESS)
      ENZO_FAIL("gFLDSplit Evolve: Error in ChemistrySource routine");

    //   solve local chemistry/gas energy systems
    if (this->AnalyticChemistry(U0, sol, extsrc, thisdt) != SUCCESS) 
      ENZO_FAIL("gFLDSplit Evolve: Error in AnalyticChemistry routine");

    // update chemistry time step size based on changes to chem+energy
//     dtchem = (this->ComputeTimeStep(U0,sol,1))*TimeUnits;
    //   (limit growth at each cycle)
    dtchem2 = this->ComputeTimeStep(U0,sol,1)*TimeUnits;
    dtchem = min(dtchem2, 2.0*dtchem);

    // enforce a solution floor on number densities
    epsilon2 = 0.0;    // put a hard floor of 0 on these fields
    if (Nchem > 0)
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
	sol_HI[i] = min(max(sol_HI[i],epsilon2),rho[i]*HFrac);
    if (Nchem > 1) {
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
	sol_HeI[i] = max(sol_HeI[i],epsilon2);
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
	sol_HeII[i] = max(sol_HeII[i],epsilon2);
    }

    //   Add fluid correction to fluid energy field (with floor)
    eh_tot = ThisGrid->GridData->AccessTotalEnergy();
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      eh_tot[i] = max(eh_tot[i]+sol_ec[i]*ecScale,tiny_number);
    if (DualEnergyFormalism) {
      eh_gas = ThisGrid->GridData->AccessGasEnergy();
      for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
	eh_gas[i] = max(eh_gas[i]+sol_ec[i]*ecScale,tiny_number);
    }

    //   Update Enzo chemistry arrays with new values
    if (Nchem > 0)  
      U0->copy_component(sol, 2);  // HI
    if (Nchem > 1) {
      U0->copy_component(sol, 3);  // HeI
      U0->copy_component(sol, 4);  // HeII
    }

    // break out of time-stepping loop if we've reached the end
    if (tchem >= tnew)  break;

  }


//   // output typical/maximum values
//   UTypVals[0] = sol->rmsnorm_component(0);
//   UMaxVals[0] = sol->infnorm_component(0);
//   for (ns=1; ns<=Nchem; ns++) {
//     UTypVals[1+ns] = U0->rmsnorm_component(ns+1);
//     UMaxVals[1+ns] = U0->infnorm_component(ns+1);
//   }
//   UTypVals[1] = 0.0;  UMaxVals[1] = 0.0;
//   for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
//     UTypVals[1] += eh_tot[i]*eh_tot[i];
//     UMaxVals[1] = max(UMaxVals[1],eh_tot[i]*eh_tot[i]);
//   }
//   UTypVals[1] = sqrt(UTypVals[1]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
//   UMaxVals[1] = sqrt(UMaxVals[1]);
// #ifdef USE_MPI
//   MPI_Allreduce(&(UMaxVals[1]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
//   UMaxVals[1] = dtmp;
//   MPI_Allreduce(&(UTypVals[1]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
//   UTypVals[1] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
// #endif
//   UTypVals[1] /= ecScale;
//   UMaxVals[1] /= ecScale;
//   if (debug) {
//     printf("   resulting internal (physical) quantities:\n");
//     printf("      Eg rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
// 	   UTypVals[0],UTypVals[0]*ErUnits0, UMaxVals[0],UMaxVals[0]*ErUnits0);
//     printf("      ec rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
// 	   UTypVals[1],UTypVals[1]*ecUnits, UMaxVals[1],UMaxVals[1]*ecUnits);
//     if (Nchem == 1) {
//       printf("     nHI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
// 	     UTypVals[2],UTypVals[2]*NiUnits0, UMaxVals[2],UMaxVals[2]*NiUnits0);
//     }
//     if (Nchem == 3) {
//       printf("    nHeI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
// 	     UTypVals[3],UTypVals[3]*NiUnits0, UMaxVals[3],UMaxVals[3]*NiUnits0);
//       printf("   nHeII rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
// 	     UTypVals[4],UTypVals[4]*NiUnits0, UMaxVals[4],UMaxVals[4]*NiUnits0);
//     }
//   }


  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // update the radiation time step size for next time step
  float RadDt = this->ComputeTimeStep(U0,sol,0);
  if (debug)  printf("   gFLDSplit time step = %g\n",RadDt);
  if (RadDt != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(RadDt);

  // update Enzo data with new values
  U0->copy_component(sol, 0);

  // scale back to Enzo units
  U0->scale_component(0,ErScale);
  U0->scale_component(1,ecScale);
  for (i=1; i<=Nchem; i++)  U0->scale_component(i+1, NiScale);

  //   Update dependent chemical species densities (ne, nHII, nHeIII) 
  //   using computed values
  if (Nchem == 1) {   // update ne, HII
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
      ne[i] = nHII[i];
    }
  }
  else if (Nchem == 3) {   // update ne, HII, HeIII
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
      nHII[i] = max(rho[i]*HFrac - nHI[i], 0.0);
      nHeIII[i] = max(rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i], 0.0);
      ne[i] = nHII[i] + nHeII[i]/4.0 + nHeIII[i]/2.0;
    }
  }

  // rescale dt, told, tnew, adot to normalized values
  dt /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;
  adot *= TimeUnits;


  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative wall time = %g\n\n",RTtime);

  // Return
  return SUCCESS;
 
}

#endif   // TRANSFER
