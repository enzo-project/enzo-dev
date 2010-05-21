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
/  Inexact Newton Forcing Term
/
/  written by: Daniel R. Reynolds
/  date:       March, 2006
/  modified1:
/
/  PURPOSE: Implements the inexact Newton forcing method given by 
/           INChoice to determine the I.N. forcing parameter.  
/           Currently, only the constant forcing term is implemented
/           (i.e. superlinear Newton convergence, with superlinearity 
/           given by forcing constant IN1)
/
************************************************************************/
#include "performance.h"
#include "InexactNewton.h"


int InexactNewtonSolver::InexactNewtonForce()
{
//   if (debug)
//     printf("Entering InexactNewtonSolver::InexactNewtonForce routine\n");

  // Compute forcing term based on inexact Newton method choice
  switch (INChoice) {
  case 1:    // choice not implemented yet
  case 2:    // choice not implemented yet
  default:   // constant forcing term (floor at 1e-12, ceil at 0.1)
//     INforce = IN1*fnorm + 1.0e-12;
//     INforce = (INforce > 0.1) ? 0.1 : INforce;
    INforce = IN1*fnorm;
//     fprintf(stdout,"   INforce: IN1 = %g, fnorm = %g, INforce = %g\n",IN1,fnorm,INforce);
    break;
  }
  
  // return success
  return SUCCESS;
}
