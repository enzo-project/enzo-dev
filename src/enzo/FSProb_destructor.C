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
/  Destructor routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified:   
/
/  PURPOSE: Frees all memory allocated for the implicit FS problem.
/
************************************************************************/
#ifdef TRANSFER
#ifdef USE_HYPRE
#include "FSProb.h"



FSProb::~FSProb()
{

//   if (debug)  printf("Entering FSProb::destructor routine\n");

  // delete HYPRE objects
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);

  // delete EnzoVectors and other internal arrays
  int i, j;
  //   EnzoVectors require deleting the structure
  delete U0;
  delete extsrc;
  delete sol;

  //   arrays require deleting the array
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(J);
  delete[] matentries;
  delete[] rhsentries;
  delete[] HYPREbuff;

  // delete boundary condition arrays
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) {
      if (BdryVals[i][j] != NULL)  delete[] BdryVals[i][j];
    }

}
#endif
#endif
