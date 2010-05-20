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
/  Gray Flux-Limited Diffusion Implicit Problem Class destructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



gFLDProblem::~gFLDProblem()
{

//   if (debug)  printf("Entering gFLDProblem::destructor routine\n");

#ifdef USE_HYPRE
  // delete HYPRE objects
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);
#endif

  // delete EnzoVectors and other internal arrays
  int i, j;
  //   InexactNewton solver requires deleting the structure
  delete INSolve;
  //   EnzoVectors require deleting the structure
  delete sol;
  delete U0;
  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete extsrc;
  delete rhs0;
  delete rhs;
  for (i=0; i<(Nchem+2); i++)
    delete L[i];

  //   arrays require deleting the array
#ifdef USE_HYPRE
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
#endif
  delete[] Ptmpvec;
  delete[] HYPREbuff;
  delete[] FluidEnergyCorrection;
  delete[] Temp;
  delete[] OpacityP;
  delete[] OpacityE;
  delete[] L;
  if ((Model >= 20) && (Model <= 29))
    delete[] MarshakParms;

  // delete boundary condition arrays
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) 
      if (EBdryVals[i][j] != NULL)  delete[] EBdryVals[i][j];

}
#endif
