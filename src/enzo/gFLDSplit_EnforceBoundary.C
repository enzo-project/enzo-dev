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
/  EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Enforces boundary conditions on a FLD problem vector.
/
/           Note: Neumann values are enforced on the first 
/                 layer of ghost zones using a first-order central 
/                 difference approximation to the first (outward-normal) 
/                 derivative.
/           Note: Since the internal radiation variables are comoving 
/                 and normalized, we renormalize the boundary conditions 
/                 as they are enforced to match the internal units.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"



int gFLDSplit::EnforceBoundary(EnzoVector *u)
{
//   if (debug)
//     printf("Entering gFLDSplit::EnforceBoundary routine\n");

  // get information about the vector u, and check against BC dims
  int i, i2, j, j2, k, k2, idx, idx2, idxbc;
  int udims[4], ugh[3][2];
  u->size(&udims[0], &udims[1], &udims[2], &udims[3], 
	  &ugh[0][0], &ugh[0][1], &ugh[1][0], 
	  &ugh[1][1], &ugh[2][0], &ugh[2][1]);
  if (udims[0] != LocDims[0]) 
    ENZO_FAIL(" EnforceBC: mismatched x0 dims");
  if (udims[1] != LocDims[1]) 
    ENZO_FAIL(" EnforceBC: mismatched x1 dims");
  if (udims[2] != LocDims[2]) 
    ENZO_FAIL(" EnforceBC: mismatched x2 dims");
  if (udims[3] != (2+Nchem)) 
    ENZO_FAIL(" EnforceBC: mismatched nspecies dims");

  // set some shortcuts for the EnzoVector dimensions, scaling
  int x0len = udims[0] + ugh[0][0] + ugh[0][1];
  int x1len = udims[1] + ugh[1][0] + ugh[1][1];

  float *udata = u->GetData(0);
  float dxa = dx[0]*LenUnits/a;
  // x0 left boundary
  //   Dirichlet
  if (OnBdry[0][0] && (BdryType[0][0]==1)) {
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<ugh[0][0]; i++) {
	  idxbc = k*LocDims[1] + j;
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	  udata[idx] = BdryVals[0][0][idxbc]/ErUnits;
	}
  }
  //   Neumann
  if (OnBdry[0][0] && (BdryType[0][0]==2)) {
    i = -1;  i2 = i+1;
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++) {
	idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	idxbc = k*LocDims[1] + j;
	udata[idx] = udata[idx2] + dxa*BdryVals[0][0][idxbc]/ErUnits;
      }
  }

  // x0 right boundary
  //   Dirichlet
  if (OnBdry[0][1] && (BdryType[0][1]==1)) {
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++)
	for (i=ArrDims[0]-ugh[0][1]; i<ArrDims[0]; i++) {
	  idxbc = k*LocDims[1] + j;
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i;
	  udata[idx] = BdryVals[0][1][idxbc]/ErUnits;
	}
  }
  //   Neumann
  if (OnBdry[0][1] && (BdryType[0][1]==2)) {
    i = LocDims[0];  i2 = i-1;
    for (k=0; k<LocDims[2]; k++)
      for (j=0; j<LocDims[1]; j++) {
	idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	idx2 = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i2+ugh[0][0];
	idxbc = k*LocDims[1] + j;
	udata[idx] = udata[idx2] + dxa*BdryVals[0][1][idxbc]/ErUnits;
      }
  }

  if (rank > 1) {
    float dya = dx[1]*LenUnits/a;
    // x1 left boundary
    //   Dirichlet
    if (OnBdry[1][0] && (BdryType[1][0]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<ugh[1][0]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = BdryVals[1][0][idxbc]/ErUnits;
	  }
    }
    //   Neumann
    if (OnBdry[1][0] && (BdryType[1][0]==2)) {
      j = -1;  j2 = j+1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dya*BdryVals[1][0][idxbc]/ErUnits;
	}
    }
      
    // x1 right boundary
    //   Dirichlet
    if (OnBdry[1][1] && (BdryType[1][1]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=ArrDims[1]-ugh[1][1]; j<ArrDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+ugh[2][0])*x1len + j)*x0len + i+ugh[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = BdryVals[1][1][idxbc]/ErUnits;
	  }
    }
    //   Neumann
    if (OnBdry[1][1] && (BdryType[1][1]==2)) {
      j = LocDims[1];  j2 = j-1;
      for (k=0; k<LocDims[2]; k++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k+ugh[2][0])*x1len + j2+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = i*LocDims[2] + k;
	  udata[idx] = udata[idx2] + dya*BdryVals[1][1][idxbc]/ErUnits;
	}
    }
  }  // end if rank > 1
     
  if (rank > 2) {
    float dza = dx[2]*LenUnits/a;
    // x2 left boundary
    //   Dirichlet
    if (OnBdry[2][0] && (BdryType[2][0]==1)) {
      for (k=0; k<ugh[2][0]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = BdryVals[2][0][idxbc]/ErUnits;
	  }
    }
    //   Neumann
    if (OnBdry[2][0] && (BdryType[2][0]==2)) {
      k = -1;  k2 = k+1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dza*BdryVals[2][0][idxbc]/ErUnits;
	}
    }
      
    // x2 right boundary
    //   Dirichlet
    if (OnBdry[2][1] && (BdryType[2][1]==1)) {
      for (k=ArrDims[2]-ugh[2][1]; k<ArrDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = (k*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = BdryVals[2][1][idxbc]/ErUnits;
	  }
    }
    //   Neumann
    if (OnBdry[2][1] && (BdryType[2][1]==2)) {
      k = LocDims[2];  k2 = k-1;
      for (j=0; j<LocDims[1]; j++)
	for (i=0; i<LocDims[0]; i++) {
	  idx = ((k+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idx2 = ((k2+ugh[2][0])*x1len + j+ugh[1][0])*x0len + i+ugh[0][0];
	  idxbc = j*LocDims[0] + i;
	  udata[idx] = udata[idx2] + dza*BdryVals[2][1][idxbc]/ErUnits;
	}
    }
  }  // end if rank > 2
      
//   if (debug)
//     printf("Exiting gFLDSplit::EnforceBoundary routine\n");

  // return success
  return SUCCESS;

}
#endif
