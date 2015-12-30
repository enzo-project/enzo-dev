/***********************************************************************
/
/  EVOLVE ROUTINE FOR COUPLED RADIATION-HYDRODYNAMICS-CHEMICAL 
/  KINETICS SYSTEM
/
/  written by: Daniel Reynolds
/  date:       November, 2006
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
 
#include "gFLDProblem.h"
#include "InexactNewton.h"
#include "CosmologyParameters.h"


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);




int gFLDProblem::Evolve(HierarchyEntry *ThisGrid, float deltat)
{

//   if (debug)  printf("Entering gFLDProblem::Evolve routine\n");

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
    ENZO_FAIL("Error in gFLDProblem_Evolve");
  }
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
    ENZO_FAIL("gFLDProblem Evolve: could not obtain velocity1");
  vy = ThisGrid->GridData->AccessVelocity2();
  if (vy == NULL) 
    ENZO_FAIL("gFLDProblem Evolve: could not obtain velocity2");
  vz = ThisGrid->GridData->AccessVelocity3();
  if (vz == NULL) 
    ENZO_FAIL("gFLDProblem Evolve: could not obtain velocity3");
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("gFLDProblem Evolve: could not obtain density");
  if (DualEnergyFormalism) {
    eh = ThisGrid->GridData->AccessGasEnergy();
    if (eh == NULL) 
      ENZO_FAIL("gFLDProblem Evolve: could not obtain fluid energy");
  }
  else {
    eh = ThisGrid->GridData->AccessTotalEnergy();
    if (eh == NULL) 
      ENZO_FAIL("gFLDProblem Evolve: could not obtain fluid energy");
  }
  RadiationEnergy = ThisGrid->GridData->AccessRadiationFrequency0();
  if (RadiationEnergy == NULL) 
    ENZO_FAIL("gFLDProblem Evolve: could not obtain Radiation energy");
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
  if ((Nchem > 0) && (nHI == NULL))
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

  // initialize external sources to 0
  extsrc->constant(0.0);

#ifdef EMISSIVITY
  // if using external Emissivity field source, copy into extsrc
  if (StarMakerEmissivityField > 0) {
    // access external emissivity field 
    float *EmissivitySource = ThisGrid->GridData->AccessEmissivity0();
    if (EmissivitySource == NULL) 
      ENZO_FAIL("gFLDProblem Evolve: could not access emissivity field");
    // access radiation source array 
    float *RadSrc = extsrc->GetData(0);
    if (RadSrc == NULL) 
      ENZO_FAIL("gFLDProblem Evolve: could not access Radiaton source");
    // copy data
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      RadSrc[i] = EmissivitySource[i];

    float srcNorm = extsrc->rmsnorm_component(0);
    if (debug)
      printf("gFLDProblem_Evolve: ||Emissivity||_rms = %g\n",srcNorm);
  }
#endif


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
      nHeIII[i] = max(0.0, rho[i]*(1.0-HFrac) - nHeI[i] - nHeII[i]);
    }
  }

  // get internal Enzo units (old time step)
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  float RadUnits;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, told) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) == FAIL) 
    ENZO_FAIL("Error in RadiationGetUnits.");
  // incorporate Enzo units with implicit solver unit scaling
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits*ecScale;
  NiUnits = DenUnits/mp*NiScale;

  // set a, adot, unit scalings to correct time-level values
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(told, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
    aUnits = 1.0/(1.0 + InitialRedshift);
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
    ENZO_FAIL("gFLDProblem Evolve: vector exchange_start error");

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
    printf("\n  gFLDProblem Evolve: current internal (physical) quantities:\n");
    printf("      Eg rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[0],UTypVals[0]*ErUnits, UMaxVals[0],UMaxVals[0]*ErUnits);
    printf("      ec rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[1],UTypVals[1]*ecUnits, UMaxVals[1],UMaxVals[1]*ecUnits);
    if (Nchem == 1) {
      printf("     nHI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[2],UTypVals[2]*NiUnits, UMaxVals[2],UMaxVals[2]*NiUnits);
    }
    if (Nchem == 3) {
      printf("    nHeI rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[3],UTypVals[3]*NiUnits, UMaxVals[3],UMaxVals[3]*NiUnits);
      printf("   nHeII rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[4],UTypVals[4]*NiUnits, UMaxVals[4],UMaxVals[4]*NiUnits);
    }
  }

  // set prepared flag to true
  prepared = true;

  // rescale dt, told, tnew, adot to physical values
  dt *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;
  adot /= TimeUnits;

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() == FAIL) 
    ENZO_FAIL("gFLDProblem Evolve: vector exchange_end error");

  // update any time-dependent boundary conditions on U0
  if (this->UpdateBoundary(U0,told,0) == FAIL) 
    ENZO_FAIL("gFLDProblem Evolve error: UpdateBoundary failure");

  // enforce boundary conditions on state U0
  if (this->EnforceBoundary(U0,0) == FAIL) 
    ENZO_FAIL("gFLDProblem_Evolve error: EnforceBoundary failure");

  // set up rhs0 for the current state, at current time t0
  if (this->ComputeRHS(rhs0, told, U0) == FAIL) 
    ENZO_FAIL("gFLDProblem Evolve: Error in ComputeRHS routine");

  // rescale dt, told, tnew, adot to normalized values
  dt /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;
  adot *= TimeUnits;

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
  NiUnits = DenUnits/mp*NiScale;

  // set a, adot to correct time-level values
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");

  // rescale dt, told, tnew, adot to physical values
  dt *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;
  adot /= TimeUnits;




  ////////////////////////////////////
  // Problem Solve Phase

  //   obtain initial guess for time-evolved solution
  sol->copy(U0);
  if (this->InitialGuess(sol) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("ERROR: InitialGuess failure!!");
  }

  //   set nonlinear solver parameters
  INSolve->SetDampedNewton(newt_linesearch);
  INSolve->SetMaxIters(newt_maxit);
  INSolve->SetInexactNewton(0, newt_INconst, 1.0);
  INSolve->SetNewtonTolerance(newt_tol);
  INSolve->SetNewtonNorm(newt_norm);
  INSolve->SetMinLinesearch(newt_MinLinesearch);


  // Call nonlinear solver to compute updated time step
  if (INSolve->Solve(this,sol) == FAIL) {
    this->Dump(sol);
    ENZO_FAIL("ERROR: INSolve failure!!");
  }
  // get Newton solver diagnostics
  int NewtIts = INSolve->GetNonlinearIterations();
  float FStep = INSolve->GetLinesearchStepsize();
  float FRes = INSolve->GetNonlinearResidual();
  if (FRes > newt_tol) {
    this->Dump(sol);
    ENZO_FAIL("ERROR: non-convergent Newton method!");
  }


  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // update the radiation time step size for next time step
  float RadDt = this->ComputeTimeStep(U0,sol,NewtIts,FStep,FRes);
  if (RadDt != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(RadDt);


  // enforce a solution floor on radiation, number densities
  //   first determine minimum floating-point roundoff
  float epsilon=1.0;
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;   // smallest representable number
  float epsilon2=0.0;
  //   enforce minimum values on solution (not yet copied into Enzo data)
  float *sol_Er = sol->GetData(0);
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    sol_Er[i] = max(sol_Er[i],epsilon);
  if (Nchem > 0) {
    float *sol_HI = sol->GetData(2);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
      sol_HI[i] = min(max(sol_HI[i],epsilon2),rho[i]*HFrac);
  }
  if (Nchem > 1) {
    float *sol_HeI = sol->GetData(3);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
      sol_HeI[i] = max(sol_HeI[i],epsilon2);
    float *sol_HeII = sol->GetData(4);
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
      sol_HeII[i] = max(sol_HeII[i],epsilon2);
  }

  // Rescale solution arrays to get back from solver to Enzo units
  sol->scale_component(0,ErScale);
  sol->scale_component(1,ecScale);
  for (i=1; i<=Nchem; i++)  sol->scale_component(i+1, NiScale);

  // Update Enzo data with new values
  //   Radiation Energy, Chemical Species and Fluid Correction are in 
  //   sol, while Enzo pointers to these fields are in U0
  U0->copy(sol);

  //   Add fluid correction to fluid energy field (with floor)
  float *eh_tot = ThisGrid->GridData->AccessTotalEnergy();
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
    eh_tot[i] = max(eh_tot[i]+FluidEnergyCorrection[i],tiny_number);
  if (DualEnergyFormalism) {
    float *eh_gas = ThisGrid->GridData->AccessGasEnergy();
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      eh_gas[i] = max(eh_gas[i]+FluidEnergyCorrection[i],tiny_number);
  }

  //   Update dependent chemical species densities (ne, nHII, nHeIII) 
  //   using computed values
  if (Nchem == 0) {        // do nothing
  }
  else if (Nchem == 1) {   // update ne, HII
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
  if (debug)  printf("RadHydro cumulative wall time = %g\n",RTtime);

  // Return
  return SUCCESS;
 
}

#endif   // TRANSFER
