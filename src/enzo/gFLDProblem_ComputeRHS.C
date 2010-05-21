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
/  Single-Group, Multi-Species Gray Flux-Limited Diffusion Implicit 
/  Problem Class RHS calculation routine.
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  August 13, 2007, by John Hayes; implemented function calls
/              calls for DiffRHS_2D (for rank = 2) and DiffRHS_1D (for rank = 1)
/
/  PURPOSE: Takes in EnzoVector and returns right-hand side of ODEs for 
/           relevant equations.  This routine will be called repeatedly, 
/           so it should NOT allocate any memory.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"
#include "CosmologyParameters.h"




int gFLDProblem::ComputeRHS(EnzoVector *rhsvec, float time, EnzoVector *u) 
{
//   if (debug)  printf("Entering gFLDProblem::ComputeRHS routine\n");

  // get local mesh description and check input vector sizes
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("ComputeRHS error: x0 vector dims do not match");
  if (usz[1] != LocDims[1])
    ENZO_FAIL("ComputeRHS error: x1 vector dims do not match");
  if (usz[2] != LocDims[2])
    ENZO_FAIL("ComputeRHS error: x2 vector dims do not match");
  if (usz[3] != (2+Nchem))
    ENZO_FAIL("ComputeRHS error: nspecies dims do not match");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("ComputeRHS error: x0 vector sizes do not match");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("ComputeRHS error: x1 vector sizes do not match");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("ComputeRHS error: x2 vector sizes do not match");

  // local rhs for all equations
  if (this->LocRHS(rhsvec, time, u) == FAIL)
    ENZO_FAIL("ComputeRHS error: LocRHS failure");

  // diffusive rhs for radiation energy equation
  tmp1->constant(0.0);   // temporary EnzoVector stored in gFLDProblem object
  float *tmp1_E = tmp1->GetData(0);
  float *Er = u->GetData(0);
  float *Er0 = U0->GetData(0);
  // note: OpacityE and Temp have already been filled in by LocRHS
  //
  if (this->DiffRHS(tmp1_E, Er, Er0, Temp, OpacityE) != SUCCESS)
    ENZO_FAIL("ComputeRHS: Error in DiffRHS routine");

  //   combine pieces together and delete temporary storage
  rhsvec->axpy_component(1.0,tmp1,0);

  // return success
  return SUCCESS;
}

#endif
