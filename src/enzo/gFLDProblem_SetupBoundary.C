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
/  Gray Flux-Limited Diffusion Implicit Problem Class SetupBoundary 
/  Routine
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  
/
/  PURPOSE: This routine is used to set values at boundary faces 
/           to be used within enforcing boundary conditions on the 
/           radiation energy field on any dimension {0,1,2}, and 
/           for either low or high face {0,1}.
/
/           For setting the face values, either scalar or array-valued 
/           arguments may be given to set the values on that face, 
/           depending on the value of BdryConst (zero implies 
/           array-valued, nonzero implies constant).
/
/           If a scalar is desired, it will be taken as the value 
/           referenced by the pointer BdryValue (pass the constant 
/           by reference).
/
/           If array input is chosen, the BdryValue array will 
/           be copied in column-major (Fortran-style) ordering:
/
/                   dim  |  fast index  |  slow index
/                 -------------------------------------
/                    0   |      x1      |      x2
/                    1   |      x2      |      x0
/                    2   |      x0      |      x1
/                 -------------------------------------
/
/           These boundary values should be given over the full extent 
/           of the domain boundary.
/
/           For Dirichlet conditions, i.e.
/               u(x) = g(x)  on boundary,
/           the value of g(x) must be supplied on cell exterior to the face.
/
/           For Neumann conditions, i.e.
/               n*grad(u(x)) = g(x)  on boundary, 
/           the flux value g(x) must be supplied on the face.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



int gFLDProblem::SetupBoundary(int Dim, int Face, 
			       int BdryConst, float *BdryData) 
{
//   if (debug)
//     printf("Entering gFLDProblem::SetupBoundary routine\n");

  // local variables
  int size, facesize, index, i;

  // Error check
  if ((Dim < 0) || (Dim >= rank)) {
    fprintf(stderr, "SetupBoundary: Dim %"ISYM" out of bounds.\n", Dim);
    ENZO_FAIL("Error in gFLDProblem_SetupBoundar");
  }
  if ((Face != 0) && (Face != 1)) {
    fprintf(stderr, "SetupBoundary: Face %"ISYM" != {0,1}.\n", Face);
    ENZO_FAIL("Error in gFLDProblem_SetupBoundar");
  }

  // compute size of local mesh and relevant faces
  size = 1;
  for (i=0; i<rank; i++)  size *= LocDims[i];
  facesize = size/LocDims[Dim];

  // delete previous boundary values arrays
  if (EBdryVals[Dim][Face] != NULL)
    delete[] EBdryVals[Dim][Face];
  
  // initialize and set boundary conditions
  EBdryVals[Dim][Face] = new float[facesize];
  
  // set constant face conditions
  if (BdryConst) {
    for (index=0; index<facesize; index++) 
      EBdryVals[Dim][Face][index] = *BdryData;
  }
  
  // set spatially-varying face conditions
  else {
    for (index=0; index<facesize; index++) 
      EBdryVals[Dim][Face][index] = BdryData[index];
  }

  return SUCCESS;
}
#endif
