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
/  Inexact Newton Solver constructor
/
/  written by: Daniel R. Reynolds
/  date:       March, 2006
/  modified1:
/
/  PURPOSE:    Deletes locally-owned temporary vectors.
/
************************************************************************/
#include "performance.h"
#include "InexactNewton.h"


InexactNewtonSolver::~InexactNewtonSolver()
{
//   if (debug)
//     printf("Entering InexactNewtonSolver::destructor\n");

  // delete local vectors tmp1, dx and fvec
  delete tmp1;
  delete dx;
  delete fvec;
}
