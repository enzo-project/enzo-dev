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
/  Gray Flux-Limited Diffusion Implicit Problem Class local rhs
/  calculation routine.
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  
/
/  PURPOSE: Takes in EnzoVector and returns local time-fixed rhs for 
/           relevant equations.  This routine will be called repeatedly, 
/           so it should NOT allocate any memory.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



int gFLDProblem::LocRHS(EnzoVector *locrhs, float time, EnzoVector *u) 
{
//   if (debug)  printf("Entering gFLDProblem::LocRHS routine\n");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("LocRHS error: x0 vector dims do not match");
  if (usz[1] != LocDims[1]) 
    ENZO_FAIL("LocRHS error: x1 vector dims do not match");
  if (usz[2] != LocDims[2]) 
    ENZO_FAIL("LocRHS error: x2 vector dims do not match");
  if (usz[3] != (2+Nchem)) 
    ENZO_FAIL("LocRHS error: nspecies dims do not match");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("LocRHS error: x0 vector sizes do not match");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("LocRHS error: x1 vector sizes do not match");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("LocRHS error: x2 vector sizes do not match");

  // extract fluid energy, radiation energy and chemistry arrays
  float *Er = u->GetData(0);
  float *ec = u->GetData(1);
  float *n_HI, *n_HeI, *n_HeII;
  if (Nchem == 0) {
    n_HI   = NULL;
    n_HeI  = NULL;
    n_HeII = NULL;
  }
  else if (Nchem == 1) {
    n_HI   = u->GetData(2);
    n_HeI  = NULL;
    n_HeII = NULL;
  }
  else if (Nchem == 3) {
    n_HI   = u->GetData(2);
    n_HeI  = u->GetData(3);
    n_HeII = u->GetData(4);
  }
  else 
    ENZO_FAIL("LocRHS ERROR: only valid for Nchem = {0,1,3}");

  // compute temperature over domain
  if (this->ComputeTemperature(Temp,time,a,u) == FAIL) 
    ENZO_FAIL("LocRHS: Error in ComputeTemperature routine");

  int size=usz[0]*usz[1]*usz[2];
  int ix, iy, iz, idx, dim;
  idx = (ghZl*ArrDims[1] + ghYl)*ArrDims[0] + ghXl;
  float *tmpArr, meanVal, minVal, maxVal;

  // compute opacity over domain
  if (this->Opacity(OpacityP, OpacityE, &time, n_HI, n_HeI, 
		    n_HeII, Temp) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in Opacity routine");

  // compute local rhs source terms over domain
  float *src_E = extsrc->GetData(0);
  float *src_e = extsrc->GetData(1);
  float *src_HI, *src_HeI, *src_HeII;
  if (Nchem == 0) {
    src_HI   = NULL;
    src_HeI  = NULL;
    src_HeII = NULL;
  }
  else if (Nchem == 1) {
    src_HI   = extsrc->GetData(2);          // Hydrogen-only case
    src_HeI  = NULL;
    src_HeII = NULL;
  }
  else if (Nchem == 3) {
    src_HI   = extsrc->GetData(2);          // Hydrogen and Helium case
    src_HeI  = extsrc->GetData(3);
    src_HeII = extsrc->GetData(4);
  }

  //    compute point-source emissivities (if not computed externally)
#ifdef EMISSIVITY
  if (StarMakerEmissivityField == 0) {
#endif
  if (this->RadiationSource(src_E, &time, Er, ec, n_HI, 
			    n_HeI, n_HeII, Temp) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in EmissivitySource routine");
#ifdef EMISSIVITY
  }
#endif

  //    compute external gas energy sources over domain
  if (this->GasEnergySource(src_e, &time, Er, ec, n_HI, 
			    n_HeI, n_HeII, Temp) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in GasEnergySource routine");

  //    compute external chemistry sources over domain
  if (this->ChemistrySource(src_HI, src_HeI, src_HeII, &time, Er, 
			    ec, n_HI, n_HeI, n_HeII, Temp) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in ChemistrySource routine");

  //    set pointers to the rhs vectors
  float *rhs_E = locrhs->GetData(0);
  float *rhs_e = locrhs->GetData(1);
  float *rhs_HI, *rhs_HeI, *rhs_HeII;
  if (Nchem == 0) {
    rhs_HI   = NULL;                        // no chemistry case
    rhs_HeI  = NULL;
    rhs_HeII = NULL;
  }
  else if (Nchem == 1) {
    rhs_HI   = locrhs->GetData(2);          // Hydrogen-only case
    rhs_HeI  = NULL;
    rhs_HeII = NULL;
  }
  else if (Nchem == 3) {
    rhs_HI   = locrhs->GetData(2);          // Hydrogen and Helium case
    rhs_HeI  = locrhs->GetData(3);
    rhs_HeII = locrhs->GetData(4);
  }

  // compute local rhs for various equations
  if (this->LocalRHS(rhs_E, rhs_e, rhs_HI, rhs_HeI, rhs_HeII, 
		     src_E, src_e, src_HI, src_HeI, src_HeII, 
		     &time, ec, Er, Temp, OpacityP, OpacityE, 
		     n_HI, n_HeI, n_HeII) != SUCCESS) 
    ENZO_FAIL("LocRHS: Error in LocalRHS routine");

  // return success
  return SUCCESS;
}

#endif
