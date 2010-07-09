/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SET A GRID'S BOUNDARY)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
// This is used to set the corners (which are not really used) of the
//   grid to something reasonable in the case of periodic B.C.'s
 
#define USE_PERIODIC
 
//
// Given a pointer to a field and its field type, find the equivalent
//   field type in the list of boundary's and apply that boundary value/type.
//   Returns: 0 on failure
//
int ExternalBoundary::SetExternalBoundary(int FieldRank, int GridDims[],
					  int GridOffset[],
					  int StartIndex[], int EndIndex[],
					  float *Field, int FieldType)
{
 
  /* declarations */
 
  int i, j, k, dim, Sign, bindex;
  float *index;
 
  /* error check: grid ranks */
 
  if (FieldRank != BoundaryRank) {
    ENZO_VFAIL("FieldRank(%"ISYM") != BoundaryRank(%"ISYM").\n",
            FieldRank, BoundaryRank)
  }
 
  /* find requested field type */
 
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (FieldType == BoundaryFieldType[field]) break;
  if (field == NumberOfBaryonFields) {
    ENZO_VFAIL("Field type (%"ISYM") not found in Boundary.\n", FieldType)
  }
 
  /* error check: make sure the boundary type array exists */
 
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
      if (BoundaryType[field][dim][0] == NULL) {
	ENZO_FAIL("BoundaryType not yet declared.\n");
      }
    }
 
  

  /* set Boundary conditions */
 
  Sign = 1;
  if (FieldType == Velocity1) Sign = -1;
 
  if (BoundaryDimension[0] > 1 && GridOffset[0] == 0) {
 
    /* set x inner (left) face */
 
    for (i = 0; i < StartIndex[0]; i++)
      for (j = 0; j < GridDims[1]; j++)
	for (k = 0; k < GridDims[2]; k++) {
	  index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
	  bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];
	  switch (BoundaryType[field][0][0][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[0] - 1 - 2*i)));
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[0]     -   i)) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][0][0][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EndIndex[0] - StartIndex[0] + 1));
#endif /* USE_PERIODIC */
	    break;
	  case shearing:
	    //	    *index = *(index + (EndIndex[0] - StartIndex[0] + 1));
	    break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (x-left).\n");
	  }
	}
  }
 
  if (BoundaryDimension[0] > 1 && GridOffset[0]+GridDims[0] == BoundaryDimension[0]) {
 
    /* set x outer (right) face */
 
    for (i = 0; i < GridDims[0]-EndIndex[0]-1; i++)
      for (j = 0; j < GridDims[1]; j++)
	for (k = 0; k < GridDims[2]; k++) {
	  index = Field + i + EndIndex[0]+1 +
	    j*GridDims[0] + k*GridDims[1]*GridDims[0];
	  bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];
	  switch (BoundaryType[field][0][1][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*i + 1)));
	    break;
	  case outflow:
	    *index =       *(index + (-1 - i)) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][0][1][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EndIndex[0] - StartIndex[0] + 1));
#endif /* USE_PERIODIC */
	    break;
	  case shearing:
	    //	    *index = *(index - (EndIndex[0] - StartIndex[0] + 1));
	    break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (x-right).\n");
	  }
	}							
  }
 
  /* set y inner (left) face */
 
  Sign = 1;
  if (FieldType == Velocity2) Sign = -1;
 
  if (BoundaryDimension[1] > 1 && GridOffset[1] == 0) {
 
    for (j = 0; j < StartIndex[1]; j++)
      for (i = 0; i < GridDims[0]; i++)
	for (k = 0; k < GridDims[2]; k++) {
	  index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
	  bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];
	  switch (BoundaryType[field][1][0][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[1] - 1 - 2*j)*GridDims[0]));
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[1]     - j)*GridDims[0]) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][1][0][bindex];
	     break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
#endif /* USE_PERIODIC */
	     break;
	  case shearing:
//	    *index = *(index + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
	     break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (y-left).\n");
	  }
	}
  }
 
  if (BoundaryDimension[1] > 1 && GridOffset[1]+GridDims[1] == BoundaryDimension[1]) {
 
    /* set y outer (right) face */
 
    for (j = 0; j < GridDims[1]-EndIndex[1]-1; j++)
      for (i = 0; i < GridDims[0]; i++)
	for (k = 0; k < GridDims[2]; k++) {
	  index = Field + i + (j + EndIndex[1]+1)*GridDims[0] +
	    k*GridDims[1]*GridDims[0];
	  bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];
	  switch (BoundaryType[field][1][1][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*j + 1)*GridDims[0]));
	    break;
	  case outflow:
	    *index =       *(index + (-1 - j)*GridDims[0]) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][1][1][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
#endif /* USE_PERIODIC */
	    break;
	  case shearing:
//	    *index = *(index - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
	    break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (y-right).\n");
	  }
	}							
  }
 
  /* set z inner (left) face */
 
  Sign = 1;
  if (FieldType == Velocity3) Sign = -1;
 
  if (BoundaryDimension[2] > 1 && GridOffset[2] == 0) {
 
    for (k = 0; k < StartIndex[2]; k++)
      for (i = 0; i < GridDims[0]; i++)
	for (j = 0; j < GridDims[1]; j++) {
	  index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
	  bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];
	  switch (BoundaryType[field][2][0][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index + (2*StartIndex[2]-1 - 2*k)*GridDims[0]*GridDims[1]));
	    break;
	  case outflow:
	    *index =       *(index + (  StartIndex[2]   - k)*GridDims[0]*GridDims[1]) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][2][0][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
#endif /* USE_PERIODIC */
	    break;
	  case shearing:
//	    *index = *(index + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
	    break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (z-left).\n");
	  }
	}
  }
 
  if (BoundaryDimension[2] > 1 && GridOffset[2]+GridDims[2] == BoundaryDimension[2]) {
 
    /* set z outer (right) face */
 
    for (k = 0; k < GridDims[2]-EndIndex[2]-1; k++)
      for (i = 0; i < GridDims[0]; i++)
	for (j = 0; j < GridDims[1]; j++) {
	  index = Field + i + j*GridDims[0] +
	    (k + EndIndex[2]+1)*GridDims[1]*GridDims[0];
	  bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];
	  switch (BoundaryType[field][2][1][bindex]) {
	  case reflecting:
	    *index = Sign*(*(index - (2*k + 1)*GridDims[0]*GridDims[1]));
	    break;
	  case outflow:
	    *index =       *(index + (-1 - k)*GridDims[0]*GridDims[1]) ;
	    break;
	  case inflow:
	    *index = BoundaryValue[field][2][1][bindex];
	    break;
	  case periodic:
#ifdef USE_PERIODIC
	    *index = *(index - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
#endif /* USE_PERIODIC */
	    break;
	  case shearing:
//	    *index = *(index - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
	    break;
	  case BoundaryUndefined:
	  default:
	    ENZO_FAIL("BoundaryType not recognized (z-right).\n");

	  }
	}							
  }
 
  return SUCCESS;
 
}
