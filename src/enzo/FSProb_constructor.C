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
/  Constructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified:   
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"

FSProb::FSProb()
{

//   if (debug)  printf("\nEntering FSProb::constructor routine\n");
  int dim, face;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif

  // initialize total solver time to zero
  FStime = 0.0;

  // initialize solver values to -1/NULL
  initial_guess = -1;
  mattype = -1;
  stSize = -1;
#ifdef USE_HYPRE
  grid = NULL;
  stencil = NULL;
#endif
  sol_tolerance = -1;
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
  maxdt = huge_number;
  tnew = -1.0;
  told = -1.0;
  dt = -1.0;
  dt_suggest = -1.0;
  theta = -1.0;
  LimType = -1;
  U0 = NULL;
  sol = NULL;
  extsrc = NULL;
  kappa = NULL;

  // initialize problem defining data 
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;
  aUnits = 1.0;
  EScale = 1.0;
  EUnits = 1.0;
  EUnits0 = 1.0;
  LenUnits = 1.0;
  LenUnits0 = 1.0;
  TimeUnits = 1.0;
  TimeUnits0 = 1.0;
  DenUnits = 1.0;
  DenUnits0 = 1.0;
  kappa0 = -1.0;
  kappa_h2on = 0;

  // initialize ionization parameters
  NGammaDot = -1.0;
  EtaRadius = -1.0;
  EtaCenter[0] = -1.0;
  EtaCenter[1] = -1.0;
  EtaCenter[2] = -1.0;

  // initialize linear solver temporary arrays to NULL
  matentries = NULL;
  rhsentries = NULL;
  HYPREbuff = NULL;

}
#endif
