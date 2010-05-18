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
/  Gray Flux-Limited Diffusion Implicit Problem linear Newton system 
/  setup function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: Called by implicit solver to notify the Problem of 
/           updates to the current state (given in the vector u), 
/           so that the linear system matrix J(u) may be updated 
/           if necessary.  For the gray FLD problem, we here compute 
/           only the local Jacobian components over the domain, and 
/           leave the actual matrix setup/solve for the lsolve 
/           routine, since we use a Schur-complement formulation for 
/           the linear system solution.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



int gFLDProblem::lsetup(EnzoVector *u)
{
//   if (debug)  printf("Entering gFLDProblem::lsetup routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) 
    ENZO_FAIL("lsetup error: gFLDProblem not yet prepared");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("lsetup error: x0 vector dims do not match");
  if (usz[1] != LocDims[1]) 
    ENZO_FAIL("lsetup error: x1 vector dims do not match");
  if (usz[2] != LocDims[2]) 
    ENZO_FAIL("lsetup error: x2 vector dims do not match");
  if (usz[3] != (2+Nchem)) 
    ENZO_FAIL("lsetup error: nspecies do not match");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("lsetup error: x0 vector sizes do not match");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("lsetup error: x1 vector sizes do not match");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("lsetup error: x2 vector sizes do not match");

  // clear Jacobian data arrays
  int i, j, k, ns, outidx;
  for (ns=0; ns<(2+Nchem); ns++)
    L[ns]->constant(0.0);


  // depending on 'Model' we either approximate or analytically compute the 
  // local Jacobian components
  //    analytically compute Jacobian  
  //    ['Models' with analytic Jacobian right now]
  if ((approx_jac == 0) && 
      ((Model == 10) || (Model == 4) || (Model == 5) || (Model == 1))) {

    // access variables, Jacobian blocks; allocate temporary variables
    float *Eg, *ec, *nHI, *nHeI, *nHeII;
    float *Egjac_Eg, *Egjac_ec, *Egjac_HI, *Egjac_HeI, *Egjac_HeII;
    float *ecjac_Eg, *ecjac_ec, *ecjac_HI, *ecjac_HeI, *ecjac_HeII;
    float *HIjac_Eg, *HIjac_ec, *HIjac_HI, *HIjac_HeI, *HIjac_HeII;
    float *HeIjac_Eg, *HeIjac_ec, *HeIjac_HI, *HeIjac_HeI, *HeIjac_HeII;
    float *HeIIjac_Eg, *HeIIjac_ec, *HeIIjac_HI, *HeIIjac_HeI, *HeIIjac_HeII;
    Eg = u->GetData(0);
    ec = u->GetData(1);
    Egjac_Eg = (L[0])->GetData(0);
    Egjac_ec = (L[0])->GetData(1);
    ecjac_Eg = (L[1])->GetData(0);
    ecjac_ec = (L[1])->GetData(1);

    if (Nchem == 0) {
      nHI = nHeI = nHeII = NULL;
      Egjac_HI = Egjac_HeI = Egjac_HeII = NULL;
      ecjac_HI = ecjac_HeI = ecjac_HeII = NULL;
      HIjac_Eg = HIjac_ec = HIjac_HI = HIjac_HeI = HIjac_HeII = NULL;
      HeIjac_Eg = HeIjac_ec = HeIjac_HI = HeIjac_HeI = HeIjac_HeII = NULL;
      HeIIjac_Eg = HeIIjac_ec = HeIIjac_HI = HeIIjac_HeI = HeIIjac_HeII = NULL;
    }
    else if (Nchem == 1) {
      nHI = u->GetData(2);
      nHeI = nHeII = NULL;      

      Egjac_HI = (L[0])->GetData(2);
      Egjac_HeI = Egjac_HeII = NULL;

      ecjac_HI = (L[1])->GetData(2);
      ecjac_HeI = ecjac_HeII = NULL;

      HIjac_Eg = (L[2])->GetData(0);
      HIjac_ec = (L[2])->GetData(1);
      HIjac_HI = (L[2])->GetData(2);
      HIjac_HeI = HIjac_HeII = NULL;

      HeIjac_Eg = HeIjac_ec = HeIjac_HI = HeIjac_HeI = HeIjac_HeII = NULL;
      HeIIjac_Eg = HeIIjac_ec = HeIIjac_HI = HeIIjac_HeI = HeIIjac_HeII = NULL;
    }
    else if (Nchem == 3) {
      nHI   = u->GetData(2);
      nHeI  = u->GetData(3);
      nHeII = u->GetData(4);

      Egjac_HI   = (L[0])->GetData(2);
      Egjac_HeI  = (L[0])->GetData(3);
      Egjac_HeII = (L[0])->GetData(4);

      ecjac_HI   = (L[1])->GetData(2);
      ecjac_HeI  = (L[1])->GetData(3);
      ecjac_HeII = (L[1])->GetData(4);

      HIjac_Eg   = (L[2])->GetData(0);
      HIjac_ec   = (L[2])->GetData(1);
      HIjac_HI   = (L[2])->GetData(2);
      HIjac_HeI  = (L[2])->GetData(3);
      HIjac_HeII = (L[2])->GetData(4);

      HeIjac_Eg   = (L[3])->GetData(0);
      HeIjac_ec   = (L[3])->GetData(1);
      HeIjac_HI   = (L[3])->GetData(2);
      HeIjac_HeI  = (L[3])->GetData(3);
      HeIjac_HeII = (L[3])->GetData(4);

      HeIIjac_Eg   = (L[4])->GetData(0);
      HeIIjac_ec   = (L[4])->GetData(1);
      HeIIjac_HI   = (L[4])->GetData(2);
      HeIIjac_HeI  = (L[4])->GetData(3);
      HeIIjac_HeII = (L[4])->GetData(4);
    }

    // call the analytic jacobian computation routine
    if (this->LocalJac(Egjac_Eg, Egjac_ec, Egjac_HI, Egjac_HeI, Egjac_HeII, 
		       ecjac_Eg, ecjac_ec, ecjac_HI, ecjac_HeI, ecjac_HeII, 
		       HIjac_Eg, HIjac_ec, HIjac_HI, HIjac_HeI, HIjac_HeII, 
		       HeIjac_Eg, HeIjac_ec, HeIjac_HI, HeIjac_HeI, 
		       HeIjac_HeII, HeIIjac_Eg, HeIIjac_ec, HeIIjac_HI, 
		       HeIIjac_HeI, HeIIjac_HeII, &tnew, Eg, ec, nHI, nHeI, 
		       nHeII) == FAIL) 
      ENZO_FAIL("lsetup: LocalJac failure!");

    // rescale local Jacobians by dt*theta for implicit time-stepping,
    // and add identity contribution
    for (i=0; i<2+Nchem; i++) {
      for (j=0; j<2+Nchem; j++) {
	(L[i])->scale_component(j,-dt*theta);
      }
      (L[i])->addconst_component(i,1.0);
    }

  }

  //    approximate Jacobian
  else {

    // we compute the local components of the Jacobian system via 
    // finite-differencing.  
    //   get typical values for input vectors
    float utypical[2+Nchem];
    utypical[0] = U0->rmsnorm_component(0);  // radiation energy
    for (i=2; i<Nchem+2; i++)                // chemistry 
      utypical[i] = U0->rmsnorm_component(i);

    //      override fluid energy correction typical value (since ec0=0)
    float dtmp1=0.0, dtmp2=0.0;
    for (i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)
      dtmp1 += eh[i]*eh[i];
    dtmp1 = sqrt(dtmp1/ArrDims[0]/ArrDims[1]/ArrDims[2])/ecScale;
#ifdef USE_MPI
    if (layout[0]*layout[1]*layout[2] == 1)
      utypical[1] = dtmp1;
    else {
      MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
      MPI_Arg one = 1;
      MPI_Allreduce(&dtmp1, &dtmp2, one, DataType, MPI_SUM, MPI_COMM_WORLD);
      utypical[1] = dtmp2/NumberOfProcessors;  // estimate based on equidistribution
    }
#else
    utypical[1] = dtmp1;
#endif

    //      just in case any of the species are not currently used
    for (i=0; i<Nchem+2; i++) 
      utypical[i] = (utypical[i] == 0.0) ? 1.0 : utypical[i];

    //   set up temporary arrays
    EnzoVector *fval = tmp1;
    EnzoVector *utmp = tmp2;
    EnzoVector *ftmp = tmp3;

    //   compute the local rhs at the current state, new time (fval)
    if (this->LocRHS(fval,tnew,u) == FAIL) 
      ENZO_FAIL("lsetup error: LocRHS failure");

    // determine floating-point roundoff
    float epsilon=1.0;
    while ((1.0 + epsilon*0.5) > 1.0)  epsilon*=0.5;
  

    // approximate the local jacobian components differently, 
    // depending on whether we are using AnalyticChem or not.
    int iz, iy, ix, idx, ns2;
    float sigma, *uarray, *utmparray, *farray, *ftmparray, *Lblock;
    float typfac = 0.001;
    if (AnalyticChem == 1) {

      /////////// First, adjust fval(2:3) to use analytical residual ///////////
      if (this->AnalyticResid(U0, u, fval, dt) == FAIL) 
	ENZO_FAIL("lsetup error: AnalyticResid failure");

      /////////// Second, handle the radiation (theta method) ///////////
      //   perturb each component of the current state (utmp),
      //   for each perturbation, compute the local rhs (ftmp);
      //   with these, approximate the Jacobian columns of E wrt that component (L)
      for (ns=0; ns<(2+Nchem); ns++) {
	
	// [re]set utmp to the current state
	utmp->copy(u);
	
	// perturb the appropriate component of utmp
	uarray = u->GetData(ns);
	utmparray = utmp->GetData(ns);
	for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	  for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	    for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	      idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	      // sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	      sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
	      utmparray[idx] += sigma;
	    }
	  }
	}
	
	// compute the local rhs due to this perturbation (ftmp)
	if (this->LocRHS(ftmp,tnew,utmp) == FAIL) 
	  ENZO_FAIL("lsetup error: LocRHS failure");
	
	// store the resulting Jacobian approximations to E components
	farray = fval->GetData(0);
	ftmparray = ftmp->GetData(0);
	Lblock = (L[0])->GetData(ns);
	for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	  for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	    for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	      idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	      // sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	      sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
	      Lblock[idx] = -dt*theta*(ftmparray[idx]-farray[idx])/sigma;
	    }  // ix loop
	  }  // iy loop
	}  // iz loop
      }  // ns loop

      // add identity contribution to Jac_EE
      (L[0])->addconst_component(0,1.0);


      /////////// Now, handle the remaining components ///////////
      for (ns=0; ns<(2+Nchem); ns++) {
	
	// [re]set utmp to the current state
	utmp->copy(u);
	
	// perturb the appropriate component of utmp
	uarray = u->GetData(ns);
	utmparray = utmp->GetData(ns);
	for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	  for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	    for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	      idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	      // sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	      sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
	      utmparray[idx] += sigma;
	    }
	  }
	}
	
	// compute the AnalyticResid due to this perturbation (ftmp)
	if (this->AnalyticResid(U0, utmp, ftmp, dt) == FAIL) 
	  ENZO_FAIL("lsetup error: AnalyticResid failure");
	
	// store the resulting Jacobian approximations
	for (ns2=1; ns2<(2+Nchem); ns2++) {
	  farray = fval->GetData(ns2);
	  ftmparray = ftmp->GetData(ns2);
	  Lblock = (L[ns2])->GetData(ns);
	  for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	    for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	      for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
		idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
		// sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
		sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
		Lblock[idx] = (ftmparray[idx]-farray[idx])/sigma;
	      }  // ix loop
	    }  // iy loop
	  }  // iz loop
	}  // ns2 loop
      }  // ns loop

    // If we're not using AnalyticChem, do this as usual
    }
    else {
    
      //   perturb each component of the current state (utmp),
      //   for each perturbation, compute the local rhs (ftmp);
      //   with these, approximate the Jacobian columns wrt that component (L)
      for (ns=0; ns<(2+Nchem); ns++) {
	
	// [re]set utmp to the current state
	utmp->copy(u);
	
	// perturb the appropriate component of utmp
	uarray = u->GetData(ns);
	utmparray = utmp->GetData(ns);
	for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	  for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	    for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
	      idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	      // sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
	      sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
	      utmparray[idx] += sigma;
	    }
	  }
	}
	
	// compute the local rhs due to this perturbation (ftmp)
	if (this->LocRHS(ftmp,tnew,utmp) == FAIL) 
	  ENZO_FAIL("lsetup error: LocRHS failure");
	
	// store the resulting Jacobian approximations
	for (ns2=0; ns2<(2+Nchem); ns2++) {
	  farray = fval->GetData(ns2);
	  ftmparray = ftmp->GetData(ns2);
	  Lblock = (L[ns2])->GetData(ns);
	  for (iz=ghZl; iz<ghZl+usz[2]; iz++) {
	    for (iy=ghYl; iy<ghYl+usz[1]; iy++) {
	      for (ix=ghXl; ix<ghXl+usz[0]; ix++) {
		idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
		// sigma = sqrt(epsilon)*max(fabs(uarray[idx]),1.0);
		sigma = sqrt(epsilon)*max(fabs(uarray[idx]),typfac*utypical[ns]);
		Lblock[idx] = -dt*theta*(ftmparray[idx]-farray[idx])/sigma;
	      }  // ix loop
	    }  // iy loop
	  }  // iz loop
	}  // ns2 loop
      }  // ns loop

      // add identity contribution to self-variable local matrices
      for (i=0; i<2+Nchem; i++)  (L[i])->addconst_component(i,1.0);

    }  // if (AnalyticChem)

  }  // end approximate Jacobian


  // in case we are using the AnalyticChem solver with Model 4 (isothermal), 
  // we need to put something in the J_{ec,ec} block to avoid numerical 
  // problems (even though it is decoupled from the other equations).
  if ((AnalyticChem == 1) && (Model == 4))  
    (L[1])->addconst_component(1,1.0);

  // return success
  return SUCCESS;
}
#endif
