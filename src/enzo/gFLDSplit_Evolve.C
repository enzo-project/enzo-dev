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
#include "CosmologyParameters.h"


//#define FAIL_ON_NAN
#define NO_FAIL_ON_NAN


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);



// This routine evolves the radiation and gas energy/chemistry equations in 
// an operator-split fashion, subcycling the physics in the following manner:
//     dt_chem <= dt_rad <= dt_hydro
// Prior to completion, the routine also updates the maximum time step the 
// overall Grid module can take to meet a maximum subcycling ratio of 
// radiation to hydrodynamics.
int gFLDSplit::Evolve(HierarchyEntry *ThisGrid, float dthydro)
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

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n gFLDSplit Evolve:\n");

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

  // Get general time-related information
  tnew = ThisGrid->GridData->ReturnTime();

  // initialize external sources to 0
  extsrc->constant(0.0);

  // access/fill radiation source array (if provided externally)
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

  // adjust chemical species HI, HeI, HeII to correspond to rho (since 
  // rho and ni are advected differently in hydro solver for some reason)
  if (this->ChemBounds(ThisGrid) != SUCCESS)
    ENZO_FAIL("gFLDSplit Evolve: chemistry bound enforcement error");

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
  if (U0->exchange_start() != SUCCESS) 
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
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) {
    UTypVals[1] += eh[i]*eh[i];
    UMaxVals[1] = max(UMaxVals[1],eh[i]*eh[i]);
  }
  UTypVals[1] = sqrt(UTypVals[1]/ArrDims[0]/ArrDims[1]/ArrDims[2]);
  UMaxVals[1] = sqrt(UMaxVals[1]);
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg one = 1;
  float dtmp;
  MPI_Allreduce(&(UMaxVals[1]), &dtmp, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  UMaxVals[1] = dtmp;
  MPI_Allreduce(&(UTypVals[1]), &dtmp, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  UTypVals[1] = dtmp/NumberOfProcessors;  // estimate based on equidistribution
#endif
  UTypVals[1] /= ecScale;
  UMaxVals[1] /= ecScale;

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("gFLDSplit Evolve: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("gFLDSplit Evolve: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits*ecScale;
  NiUnits = (Nchem == 0) ? NiScale : DenUnits/mp*NiScale;
  if (debug) {
    printf("   current internal (physical) quantities:\n");
    printf("      Eg rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[0],UTypVals[0]*ErUnits, UMaxVals[0],UMaxVals[0]*ErUnits);
    printf("      ec rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	   UTypVals[1],UTypVals[1]*ecUnits, UMaxVals[1],UMaxVals[1]*ecUnits);
    if (Nchem >= 1) {
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

  // if autoScale enabled, determine scaling factor updates here
  float ScaleCorrTol = 1.e-2;
  float ErScaleCorr = 1.0;
  float ecScaleCorr = 1.0;
  float NiScaleCorr = 1.0;
  if (StartAutoScale && autoScale) {
    if ((UMaxVals[0] - UTypVals[0]) > ScaleCorrTol*UMaxVals[0])
      ErScaleCorr = UMaxVals[0];
    if ((UMaxVals[1] - UTypVals[1]) > ScaleCorrTol*UMaxVals[1])
      ecScaleCorr = UMaxVals[1];
    if (Nchem > 0) {
      float NiMax = UMaxVals[2];
      float NiTyp = UTypVals[2];
      if (Nchem > 1) {
        NiMax = max(NiMax, UMaxVals[3]);
        NiTyp = max(NiTyp, UTypVals[3]);
        NiMax = max(NiMax, UMaxVals[4]);
        NiTyp = max(NiTyp, UTypVals[4]);
      }
      if ((NiMax - NiTyp) > ScaleCorrTol*NiMax)
       NiScaleCorr = NiMax;
    }
  }

  // initialize variables that we'll use throughout the time subcycling
  float stime2, ftime2, stime3, ftime3;   // radiation, chemistry timers
  float thisdt, dtchem2, tchem;   // chemistry time-stepping variables
  int radstep, chemstep, radstop, chemstop;  // subcycle iterators
  float *sol_ec = sol->GetData(1);

  // have U0 finish communication of neighbor information
  if (U0->exchange_end() != SUCCESS) 
    ENZO_FAIL("gFLDSplit Evolve: vector exchange_end error");



  ////////////////////////////////////
  // Problem Solve Phase

  // internal time-stepping loop to catch up with Hydro time
  float end_time = tnew + dthydro;
  radstop = 0;
  for (radstep=0; radstep<=maxsubcycles*100; radstep++) {
      
    // start MPI timer for radiation solver
#ifdef USE_MPI
    stime2 = MPI_Wtime();
#else
    stime2 = 0.0;
#endif

    // update time-step information
    told = tnew;

    // keep trying time steps until radiation solver succeeds. 
    // Note: if we reach the minimum time step size, RadStep will call ENZO_FAIL
    int recompute_step = 1;
    while (recompute_step) {

      // update time-step information.  Note: dtrad was set on previous 
      // iteration of solver, or by user input for first iteration
      tnew = told + dtrad;
      if ((tnew - end_time)/end_time > -1.0e-14) {   // do not exceed synchronization time
	tnew = end_time;
	radstop = 1;
      }
      dt = tnew - told;
      if (debug) 
       printf("\n subcycled rad %"ISYM": dt=%7.1e, t=%7.1e (hydro dt=%7.1e, t=%7.1e)\n",
        radstep,dt,tnew,dthydro,end_time);
      
      // take a radiation step
      recompute_step = this->RadStep(ThisGrid, eta_set);
      
      // if the radiation step was unsuccessful, back-track to previous 
      // step and pull back on dtrad
      if (recompute_step) {
       dtrad = max(dtrad*0.5, mindt);
       tnew = told;
       radstop = 0;
      }

    }
    
    
    // stop MPI timer for radiation solver, increment total
#ifdef USE_MPI
    ftime2 = MPI_Wtime();
#else
    ftime2 = 0.0;
#endif
    HYPREtime += ftime2-stime2;
    
    
    // subcycle the chemistry and gas energy equations (if not handled elsewhere)
    if (!RadiativeCooling) {
      
      // start MPI timer for chemistry/gas solver
#ifdef USE_MPI
      stime3 = MPI_Wtime();
#else
      stime3 = 0.0;
#endif
      
      tchem = told;
      chemstop = 0;
      for (chemstep=0; chemstep<=maxchemsub*2; chemstep++) {
      	// update tchem
        thisdt = min(dtchem, dt);             // do not exceed radiation dt
        thisdt = max(thisdt, dt/maxchemsub);  // set max subcycle count wrt radiation
        tchem += thisdt;                      // update chemistry time
        if ((tchem - tnew)/tnew > -1.0e-14) { // do not exceed radiation time
          thisdt = tnew - (tchem - thisdt);   // get max time step
          tchem = tnew;                       // set updated time
          chemstop = 1;
        }
        if (debug) 
          printf("   subcycled chem %"ISYM": dt=%7.1e, t=%7.1e (rad dt=%7.1e, t=%7.1e)\n",
            chemstep,thisdt,tchem,dt,tnew);

        //   take a chemistry step
        if (this->ChemStep(ThisGrid, thisdt, tchem) != SUCCESS)
          ENZO_FAIL("gFLDSplit Evolve: Error in ChemStep routine");

        // update chemistry time step size based on changes to chem+energy
        //   (limit growth at each cycle)
        dtchem2 = this->ComputeTimeStep(U0,sol,1);
        dtchem = min(dtchem2, 2.0*dtchem);

        //   Zero out fluid energy correction fields
        for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
          sol_ec[i] = 0.0;
        for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
          FluidEnergyCorrection[i] = 0.0;

        //   Update Enzo chemistry arrays with new values
        if (Nchem > 0)  
          U0->copy_component(sol, 2);  // HI
        if (Nchem > 1) {
	        U0->copy_component(sol, 3);  // HeI
	        U0->copy_component(sol, 4);  // HeII
        }

        // break out of time-stepping loop if we've reached the end
        if (chemstop)  break;
	
      }  // end of chemistry subcycling
      
      // stop MPI timer for chemistry/gas solver, increment total
#ifdef USE_MPI
      ftime3 = MPI_Wtime();
#else
      ftime3 = 0.0;
#endif
      ChemTime += ftime3-stime3;
      
    }  // if (!RadiativeCooling)
    
    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est = this->ComputeTimeStep(U0,sol,0);
    dtrad = min(dt_est, dtgrowth*dtrad);

    // update Enzo radiation field with new values
    U0->copy_component(sol, 0);
    
    // have U0 communicate neighbor information
    if (U0->exchange() != SUCCESS) 
      ENZO_FAIL("gFLDSplit Evolve: vector exchange error");

    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop



  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // if chemistry and cooling handled elsewhere, fill the rates
  if (RadiativeCooling) {
    float *phHI       = ThisGrid->GridData->AccessKPhHI();
    float *phHeI      = ThisGrid->GridData->AccessKPhHeI();
    float *phHeII     = ThisGrid->GridData->AccessKPhHeII();
    float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
    float *dissH2I    = ThisGrid->GridData->AccessKDissH2I();
    this->FillRates(sol, U0, phHI, phHeI, phHeII, photogamma, dissH2I);
  }

  // update the radiation time step size for next time step
  if (dtrad != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(dtrad*maxsubcycles);

  // scale back to Enzo units
  U0->scale_component(0,ErScale);
  U0->scale_component(1,ecScale);
  for (i=1; i<=Nchem; i++)  U0->scale_component(i+1, NiScale);

  // update scaling factors to account for new values
  if (StartAutoScale && autoScale) {
    ErScale *= max(ErScaleCorr, 1.0);
    ecScale *= max(ecScaleCorr, 1.0);
    NiScale *= max(NiScaleCorr, 1.0);
  }

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

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative time = %g (HYPRE = %g,  Chem = %g)\n\n",
		     RTtime, HYPREtime, ChemTime);

  // Return
  return SUCCESS;
 
}



// This routine evolves one time step of the chemistry and gas energy subsystem 
// within the gFLDSplit module.  
int gFLDSplit::ChemStep(HierarchyEntry *ThisGrid, float thisdt, float tcur)
{

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tcur-thisdt) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_ChemStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tcur-thisdt) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_ChemStep: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits0 = RadUnits*ErScale;
  NiUnits0 = (Nchem == 0) ? NiScale : DenUnits0/mp*NiScale;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tcur-thisdt, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("gFLDSplit_ChemStep: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tcur) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_ChemStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tcur) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_ChemStep: Error in RadiationGetUnits.");
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits*ecScale;
  NiUnits = (Nchem == 0) ? NiScale : DenUnits/mp*NiScale;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tcur, &a, &adot) != SUCCESS) 
      ENZO_FAIL("gFLDSplit_ChemStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale thisdt, tcur to physical values for use within solver
  thisdt *= TimeUnits;
  tcur   *= TimeUnits;

  // set pointers to relevant data arrays
  float *ecsrc = extsrc->GetData(1);
  float *HIsrc   = NULL;
  float *HeIsrc  = NULL;
  float *HeIIsrc = NULL;
  float *sol_ec   = sol->GetData(1);
  float *sol_HI   = NULL;
  float *sol_HeI  = NULL;
  float *sol_HeII = NULL;
  float *eh_tot = ThisGrid->GridData->AccessTotalEnergy();
  float *eh_gas = NULL;
  if (DualEnergyFormalism)  
    eh_gas = ThisGrid->GridData->AccessGasEnergy();
  if (Nchem > 0) {
    sol_HI = sol->GetData(2);
    HIsrc  = extsrc->GetData(2);
  }
  if (Nchem > 1) {
    sol_HeI  = sol->GetData(3);
    sol_HeII = sol->GetData(4);
    HeIsrc   = extsrc->GetData(3);
    HeIIsrc  = extsrc->GetData(4);
  }
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL) 
    ENZO_FAIL("gFLDSplit_Chemstep: could not obtain density");

  //   fill in the gas energy and chemistry source terms
  if (this->GasEnergySource(ecsrc, &tcur) != SUCCESS)
    ENZO_FAIL("gFLDSplit_Chemstep: Error in GasEnergySource routine");
  if (this->ChemistrySource(HIsrc, HeIsrc, HeIIsrc, &tcur) != SUCCESS)
    ENZO_FAIL("gFLDSplit_Chemstep: Error in ChemistrySource routine");
  
  //   solve local chemistry/gas energy systems
  if (this->AnalyticChemistry(U0, sol, extsrc, thisdt) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_Chemstep: Error in AnalyticChemistry routine");
  
  // enforce a solution floor on number densities
  float epsilon2 = 0.0;    // put a hard floor of 0 on these fields
  int i;
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
  for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
    eh_tot[i] = max(eh_tot[i]+sol_ec[i]*ecScale,tiny_number);
  if (DualEnergyFormalism) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      eh_gas[i] = max(eh_gas[i]+sol_ec[i]*ecScale,tiny_number);
  }

  // rescale thisdt, tcur back to normalized values
  thisdt /= TimeUnits;
  tcur   /= TimeUnits;

  return SUCCESS;
}



// This routine adjusts Enzo's chemical species HI, HeI, HeII to correspond 
// to the density, since rho and ni are advected differently in hydro solver, 
// which can destroy cell-local mass conservation.
int gFLDSplit::ChemBounds(HierarchyEntry *ThisGrid)
{

  // Set pointers to each variable
  float *nHI    = NULL;
  float *nHII   = NULL;
  float *nHeI   = NULL;
  float *nHeII  = NULL;
  float *nHeIII = NULL;
  float *ne     = NULL;
  rho = ThisGrid->GridData->AccessDensity();
  if (rho == NULL)  ENZO_FAIL("gFLDSplit_ChemBounds: could not obtain density");
  nHI    = ThisGrid->GridData->AccessHIDensity();
  nHII   = ThisGrid->GridData->AccessHIIDensity();
  nHeI   = ThisGrid->GridData->AccessHeIDensity();
  nHeII  = ThisGrid->GridData->AccessHeIIDensity();
  nHeIII = ThisGrid->GridData->AccessHeIIIDensity();
  ne     = ThisGrid->GridData->AccessElectronDensity();
  //    check that we accessed the required species for Nchem
  if (Nchem > 0) 
    if (nHI == NULL) 
      ENZO_FAIL("gFLDSplit_ChemBounds error: cannot access HI density!");
  if (Nchem > 1) {
    if (nHeI == NULL) 
      ENZO_FAIL("gFLDSplit_ChemBounds error: cannot access HeI density!");
    if (nHeII == NULL) 
      ENZO_FAIL("gFLDSplit_ChemBounds error: cannot access HeII density!");
  }

  // Hydrogen chemistry 
  int i;
  float rhochem;
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

  // electron density
  if (Nchem == 1) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  ne[i] = nHII[i];
  }
  else if (Nchem == 3) {
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++) 
      ne[i] = nHII[i] + nHeII[i]/4.0 + nHeIII[i]/2.0;
  }

  return SUCCESS;
}



// This routine evolves the radiation subsystem within the gFLDSplit module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
int gFLDSplit::RadStep(HierarchyEntry *ThisGrid, int eta_set)
{
   
  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in RadiationGetUnits.");
  float mp = 1.67262171e-24;   // Mass of a proton [g]
  ErUnits0 = RadUnits*ErScale;
  NiUnits0 = (Nchem == 0) ? NiScale : DenUnits0/mp*NiScale;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("gFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in RadiationGetUnits.");
  ErUnits = RadUnits*ErScale;
  ecUnits = VelUnits*VelUnits*ecScale;
  NiUnits = (Nchem == 0) ? NiScale : DenUnits/mp*NiScale;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("gFLDSplit_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew, adot to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;

  // compute emissivity at this internal time step (if provided internally)
  float *RadSrc = extsrc->GetData(0);
  if (eta_set == 0) {
    if (this->RadiationSource(RadSrc, &tnew) != SUCCESS) 
      ENZO_FAIL("gFLDSplit_RadStep: Error in RadiationSource routine");
  }

  // compute emissivity statistics  
  float srcNorm = extsrc->rmsnorm_component(0);
  float srcMax  = extsrc->infnorm_component(0);
  if (debug)  printf("   emissivity norm = %g,  max = %g\n",srcNorm,srcMax);
    
  // turn on automatic scaling for next step if this step has nontrivial emissivity
  float ScaleCorrTol = 1.e-2;
  if ((srcMax - srcNorm) > ScaleCorrTol*srcMax)  
    StartAutoScale = true;

  //   enforce boundary conditions on current time step vector
  if (this->EnforceBoundary(U0) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: EnforceBoundary failure!!");
    
  //   set initial guess as old solution
  sol->copy(U0);
    
  //   compute updated opacities
  if (this->Opacity(OpacityE, &tnew, sol) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in Opacity routine");
    
  //   compute the gas 'Temperature' at old and new time steps
  if (this->ComputeTemperature(Temperature0, U0) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in ComputeTemperature routine");
  if (this->ComputeTemperature(Temperature, sol) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in ComputeTemperature routine");
    
#ifdef USE_HYPRE
  
  // set up and solve radiation equation
  float *RadiationEnergy = U0->GetData(0);    // old radiation energy array
  float *Eg_new = sol->GetData(0);    // updated radiation energy array
  float rhsnorm;   // used for setting HYPRE solver tolerance
  if (this->SetupSystem(matentries, rhsentries, &rhsnorm, RadiationEnergy, Eg_new, 
			OpacityE, Temperature, Temperature0, RadSrc) != SUCCESS) 
    ENZO_FAIL("gFLDSplit_RadStep: Error in SetupSystem routine");
  
  // skip solve if ||rhs|| < sol_tolerance  (i.e. old solution is fine)
  if (rhsnorm < sol_tolerance) {
    if (debug) {
      printf(" ----------------------------------------------------------------------\n");
      printf("   no solve required: |rhs| = %.1e  <  tol = %.1e\n", rhsnorm, sol_tolerance);
      printf(" ======================================================================\n\n");
    }
    // rescale dt, told, tnew, adot back to normalized values
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;
    return 0;
  }
  
  // assemble matrix
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};   // matrix stencil entries
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  Eflt64 delta;    // used for setting HYPRE solver tolerance
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, matentries); 
  HYPRE_StructMatrixAssemble(P);
  
  // insert rhs into HYPRE vector b
  HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, rhsentries);
  
  // set linear solver tolerance (rescale to relative residual and not actual)
  delta = (rhsnorm > 1.e-8) ? sol_tolerance/rhsnorm : delta;
  //  delta = min(delta, 1.0e-6);
  delta = min(delta, 1.0e-2);
  
  // insert sol initial guess into HYPRE vector x 
  int xBuff, yBuff, zBuff, Zbl, Ybl, ix, iy, iz;  // mesh indexing shortcuts
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  xBuff = GhDims[0][0]-SolvOff[0];
  yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
  zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];  ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];  ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = 0.0;
      HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
    }
  }
  
  // assemble vectors
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

  // set up the solver and preconditioner [PFMG]
  //    create the solver & preconditioner
  HYPRE_StructSolver solver;            // HYPRE solver structure
  HYPRE_StructSolver preconditioner;    // HYPRE preconditioner structure
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    break;
  }
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);
  
  // Multigrid solver: for periodic dims, only coarsen until grid no longer divisible by 2
  Eint32 max_levels, level=-1;
  int Ndir;
  if (BdryType[0][0] == 0) {
    level = 0;
    Ndir = GlobDims[0];
    while ( Ndir%2 == 0 ) {
      level++;
      Ndir /= 2;
    }
  }
  max_levels = level;
  if (rank > 1) {
    if (BdryType[1][0] == 0) {
      level = 0;
      Ndir = GlobDims[1];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  if (rank > 2) {
    if (BdryType[2][0] == 0) {
      level = 0;
      Ndir = GlobDims[2];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }

  //    set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_maxit/4);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  
  //    set solver options
  switch (Krylov_method) {
  case 0:   // PCG
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
    else {    // ignore preconditioner for 1D tests (bug); increase CG its
      HYPRE_StructPCGSetMaxIter(solver, sol_maxit*500);
    }
    if (delta != 0.0)   HYPRE_StructPCGSetTol(solver, delta);
    HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
    break;
  case 2:   // GMRES
    //  HYPRE_StructGMRESSetPrintLevel(solver, sol_printl);
    HYPRE_StructGMRESSetLogging(solver, sol_log);
    //  HYPRE_StructGMRESSetRelChange(solver, 1);
    if (rank > 1) {
      HYPRE_StructGMRESSetMaxIter(solver, sol_maxit);
      HYPRE_StructGMRESSetKDim(solver, sol_maxit);
      HYPRE_StructGMRESSetPrecond(solver, 
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				  preconditioner);
    }
    else {    // ignore preconditioner for 1D tests (bug); increase CG its
      HYPRE_StructGMRESSetMaxIter(solver, sol_maxit*50);
      HYPRE_StructGMRESSetKDim(solver, sol_maxit*50);
    }
    if (delta != 0.0)   HYPRE_StructGMRESSetTol(solver, delta);
    HYPRE_StructGMRESSetup(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    //  HYPRE_StructBiCGSTABSetPrintLevel(solver, sol_printl);
    HYPRE_StructBiCGSTABSetLogging(solver, sol_log);
    if (rank > 1) {
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_maxit);
      HYPRE_StructBiCGSTABSetPrecond(solver, 
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				     preconditioner);
    }
    else {    // ignore preconditioner for 1D tests (bug); increase its
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_maxit*500);
    }
    if (delta != 0.0)   HYPRE_StructBiCGSTABSetTol(solver, delta);
    HYPRE_StructBiCGSTABSetup(solver, P, rhsvec, solvec);
    break;
  }
  
  // solve the linear system
  if (debug)
    printf(" ----------------------------------------------------------------------\n");
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESSolve(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABSolve(solver, P, rhsvec, solvec);
    break;
  }
  
  // extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;  // HYPRE solver statistics
  Eint32 Sits=0, Pits=0;  // HYPRE solver statistics
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructPCGGetNumIterations(solver, &Sits);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructGMRESGetNumIterations(solver, &Sits);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructBiCGSTABGetNumIterations(solver, &Sits);
    break;
  }
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug) printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
		    finalresid, delta, rhsnorm, Sits, Pits);
  int recompute_step = 0;
  if ((sol_tolerance != 0.0) || (finalresid != finalresid)) {
//   if ((sol_tolerance != 0.0) || isnan(finalresid)) {
    // if final residual is too large (or nan) set return value to reduce 
    // dt and recompute step, unless we're at the minimum step size already
    if ((finalresid*rhsnorm > sol_tolerance) || (finalresid != finalresid)) {
//     if ((finalresid*rhsnorm > sol_tolerance) || isnan(finalresid)) {

#ifndef FAIL_ON_NAN
      if (dt > mindt*TimeUnits) {
	// allow remainder of function to complete (to reset units, etc.), 
	// but have calling routine update dt and compute step again.
	recompute_step = 1;
      }
      else {
#endif
	fprintf(stderr,"gFLDSplit_RadStep: could not achieve prescribed tolerance!\n");
	
	// output linear system to disk
	if (debug)  printf("Writing out matrix to file P.mat\n");
	HYPRE_StructMatrixPrint("P.mat",P,0);
	if (debug)  printf("Writing out rhs to file b.vec\n");
	HYPRE_StructVectorPrint("b.vec",rhsvec,0);
	if (debug)  printf("Writing out current solution to file x.vec\n");
	HYPRE_StructVectorPrint("x.vec",solvec,0);
	
	// dump module parameters to disk
	this->Dump(sol);
	ENZO_FAIL("Error in gFLDSplit_RadStep");
#ifndef FAIL_ON_NAN
      }
#endif
    }
  }
  if (debug)  printf(" ======================================================================\n\n");
  
  
  // if solve was successful: extract values and add to current solution
  if (!recompute_step) 
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
  
  // destroy HYPRE solver & preconditioner structures
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGDestroy(solver);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESDestroy(solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABDestroy(solver);
    break;
  }
  HYPRE_StructPFMGDestroy(preconditioner);
  
  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  for (int i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    Eg_new[i] = max(Eg_new[i],epsilon);

  // rescale dt, told, tnew, adot back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

#else
  ENZO_FAIL("gFLDSplit_RadStep ERROR: module requires USE_HYPRE to be set!");
#endif

}




#endif   // TRANSFER
