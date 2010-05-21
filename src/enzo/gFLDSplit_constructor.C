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
/  Constructor routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"


gFLDSplit::gFLDSplit()
{

//   if (debug)  printf("\nEntering gFLDSplit::constructor routine\n");
  int dim, face;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif

  // initialize total RT time to zero
  RTtime = 0.0;

  // initialize HYPRE values to -1/NULL
  mattype = -1;
  stSize = -1;
#ifdef USE_HYPRE
  grid = NULL;
  stencil = NULL;
#endif
  sol_tolerance = -1.0;
  sol_maxit = -1;
  sol_rlxtype = -1;
  sol_npre = -1;
  sol_npost = -1;
  sol_printl = -1;
  sol_log = -1;
  totIters = -1;
  for (dim=0; dim<3; dim++) {
    for (face=0; face<2; face++)
      SolvIndices[dim][face] = 0;
    SolvOff[dim] = 0;
  }

  // initialize solver values to -1
  initial_guess = -1;

  // initialize problem grid information to -1/NULL
  rank = -1;
  for (dim=0; dim<3; dim++) {
    layout[dim] = 1;    // initialize for single-processor run
    location[dim] = 0;  // initialize for single-processor run
    LocDims[dim] = 1;
    ArrDims[dim] = 1;
    dx[dim] = 1.0;
    for (face=0; face<2; face++) {
      OnBdry[dim][face] = false;
      NBors[dim][face] = MPI_PROC_NULL;
      GhDims[dim][face] = 0;
      BdryType[dim][face] = -1;
      EdgeVals[dim][face] = -1.0;
      BdryVals[dim][face] = NULL;
    }
  }
  
  // initialize time-stepping related data to -1/NULL
  maxdt = 1.0e20;
  mindt = 0.0;
  initdt = 1.0e20;
  dtfac[0] = 1.0e20;
  dtfac[1] = 1.0e20;
  dtfac[2] = 1.0e20;
  dtnorm = 0.0;
  tnew = -1.0;
  told = -1.0;
  dt = -1.0;
  dtchem = -1.0;
  theta = -1.0;
  sol = NULL;
  U0 = NULL;
  extsrc = NULL;
  

  // initialize problem defining data 
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;
  aUnits = 1.0;
  ErScale = 1.0;
  ecScale = 1.0;
  NiScale = 1.0;
  ErUnits = 1.0;
  ErUnits0 = 1.0;
  ecUnits = 1.0;
  NiUnits = 1.0;
  NiUnits0 = 1.0;
  DenUnits = 1.0;
  DenUnits0 = 1.0;
  LenUnits = 1.0;
  LenUnits0 = 1.0;
  TimeUnits = 1.0;
  VelUnits = 1.0;
  Nchem = 0;
  HFrac = 1.0;  // default is pure hydrogen (no helium)
  Model = -1;
  ESpectrum = -1;
  intSigE = 0.0;
  intSigESigHI = 0.0;
  intSigESigHeI = 0.0;
  intSigESigHeII = 0.0;
  intSigESigHInu = 0.0;
  intSigESigHeInu = 0.0;
  intSigESigHeIInu = 0.0;

  // initialize linear solver/Jacobian arrays to NULL
  matentries = NULL;
  rhsentries = NULL;
  HYPREbuff  = NULL;

  // initialize HYPRE structures to NULL
#ifdef USE_HYPRE
  P      = NULL;
  rhsvec = NULL;
  solvec = NULL;
#endif

  // initialize access to Enzo arrays to NULL
  vx = NULL;
  vy = NULL;
  vz = NULL;
  rho = NULL;
  eh = NULL;

  // initialize storage arrays to NULL
  Temperature = NULL;
  Temperature0 = NULL;
  FluidEnergyCorrection = NULL;
  OpacityE = NULL;

}
#endif
