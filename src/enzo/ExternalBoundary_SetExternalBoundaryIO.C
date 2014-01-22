/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SET A GRID'S BOUNDARY)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       November, 2005
/              Out-of-core handling for the boundary
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <hdf5.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
void my_exit(int status);

// HDF5 function prototypes



int READ_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
int READ_BV(float         *bv_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);

 
// This is used to set the corners (which are not really used) of the
//   grid to something reasonable in the case of periodic B.C.'s
 
//#define USE_PERIODIC
 
// Given a pointer to a field and its field type, find the equivalent
//   field type in the list of boundary's and apply that boundary value/type.
//   Returns: 0 on failure

int ExternalBoundary::SetExternalBoundary(int FieldRank, int GridDims[],
					  int GridOffset[],
					  int StartIndex[], int EndIndex[],
					  float *Field, int FieldType)
{
 
  /* declarations */
 
  int i, j, k, dim, Sign, bindex;
  float *index;

#ifdef OOC_BOUNDARY
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;

  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];

  hssize_t    mem_offset;
  hssize_t    file_offset[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int face;
  int slabsize;
  int cubesize;
  int MaxFaceSize;
  float *bv_buffer;
  boundary_type *bt_buffer;
  boundary_type btb;
#endif
 
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

  // for ExternalBoundaryIO it should NOT exist!

#ifndef OOC_BOUNDARY
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
      if (BoundaryType[field][dim][0] == NULL) {
	ENZO_VFAIL("BoundaryType not yet declared for field: %i.\n", field)
      }
    }
#endif
 
  /* set Boundary conditions */

  // call is by Field - set all 6 faces
  // dim = 0, face = 0
  // dim = 0, face = 1
  // dim = 1, face = 0
  // dim = 1, face = 1
  // dim = 2, face = 0
  // dim = 2, face = 1

#ifdef OOC_BOUNDARY
    cubesize = 1;
    for ( i = 0; i < BoundaryRank; i++ ) {
      cubesize = cubesize * BoundaryDimension[i];
    }

    MaxFaceSize = 1;
    for ( i = 0; i < BoundaryRank; i++ ) {
      MaxFaceSize = max(MaxFaceSize, cubesize/BoundaryDimension[i]);
    }

//    fprintf(stderr, "Assign new boundary_type bt_buffer[%"ISYM"]\n", MaxFaceSize);
    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        bt_buffer = new boundary_type[MaxFaceSize];
      }
    }

    if (ExternalBoundaryValueIO)
      bv_buffer = new float[MaxFaceSize];

#endif

  Sign = 1;
  if (FieldType == Velocity1) Sign = -1;
 
  if (BoundaryDimension[0] > 1 && GridOffset[0] == 0) {
 
    /* set x inner (left) face */

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 0;
    face = 0;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (i = 0; i < StartIndex[0]; i++)
      for (j = 0; j < GridDims[1]; j++)
        for (k = 0; k < GridDims[2]; k++) {

          index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
          bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "xl %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index + (2*StartIndex[0] - 1 - 2*i)));
            break;
          case outflow:
            *index =       *(index + (  StartIndex[0]     -   i)) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
            break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index + (EndIndex[0] - StartIndex[0] + 1));
#endif /* USE_PERIODIC */
            break;
  	  case shearing:
	     *index = *(index + (EndIndex[0] - StartIndex[0] + 1));
 	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (x-left).\n",btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType!\n");
          }

        }

#else
 
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
 	    *index = *(index + (EndIndex[0] - StartIndex[0] + 1));
	    break;
	  case BoundaryUndefined:
            break;
	  default:
	    ENZO_VFAIL("BoundaryType %"ISYM" not recognized (x-left).\n",
		    BoundaryType[field][0][0][bindex])
	  }

	}

#endif

  }
 
  if (BoundaryDimension[0] > 1 && GridOffset[0]+GridDims[0] == BoundaryDimension[0]) {
 
    /* set x outer (right) face */

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 0;
    face = 1;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (i = 0; i < GridDims[0]-EndIndex[0]-1; i++)
      for (j = 0; j < GridDims[1]; j++)
        for (k = 0; k < GridDims[2]; k++) {

          index = Field + i + EndIndex[0]+1 +
            j*GridDims[0] + k*GridDims[1]*GridDims[0];
          bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "xr %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index - (2*i + 1)));
            break;
          case outflow:
            *index =       *(index + (-1 - i)) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
            break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index - (EndIndex[0] - StartIndex[0] + 1));
#endif /* USE_PERIODIC */
            break;
 	  case shearing:
 	     *index = *(index - (EndIndex[0] - StartIndex[0] + 1));
 	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (x-right).\n",btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType!\n");
          }

        }

#else
 
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
 	    *index = *(index - (EndIndex[0] - StartIndex[0] + 1));
 	    break;
	  case BoundaryUndefined:
            break;
	  default:
	    ENZO_VFAIL("BoundaryType %"ISYM" not recognized (x-right).\n",
		    BoundaryType[field][0][1][bindex])
	  }

	}

#endif

  }
 
  /* set y inner (left) face */
 
  Sign = 1;
  if (FieldType == Velocity2) Sign = -1;
 
  if (BoundaryDimension[1] > 1 && GridOffset[1] == 0) {

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 1;
    face = 0;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (j = 0; j < StartIndex[1]; j++)
      for (i = 0; i < GridDims[0]; i++)
        for (k = 0; k < GridDims[2]; k++) {

          index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
          bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "yl %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index + (2*StartIndex[1] - 1 - 2*j)*GridDims[0]));
            break;
          case outflow:
            *index =       *(index + (  StartIndex[1]     - j)*GridDims[0]) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
             break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
#endif /* USE_PERIODIC */
             break;
 	  case shearing:
 	    *index = *(index + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (y-left).\n",
		    btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType!\n");
          }

        }

#else
 
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
 	    *index = *(index + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
 	    break;
	  case BoundaryUndefined:
            break;
	  default:
	    ENZO_VFAIL("BoundaryType %"ISYM" not recognized (y-left).\n",
		    BoundaryType[field][1][0][bindex])
	  }

	}

#endif

  }
 
  if (BoundaryDimension[1] > 1 && GridOffset[1]+GridDims[1] == BoundaryDimension[1]) {
 
    /* set y outer (right) face */

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 1;
    face = 1;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (j = 0; j < GridDims[1]-EndIndex[1]-1; j++)
      for (i = 0; i < GridDims[0]; i++)
        for (k = 0; k < GridDims[2]; k++) {

          index = Field + i + (j + EndIndex[1]+1)*GridDims[0] +
            k*GridDims[1]*GridDims[0];
          bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "yr %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index - (2*j + 1)*GridDims[0]));
            break;
          case outflow:
            *index =       *(index + (-1 - j)*GridDims[0]) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
            break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
#endif /* USE_PERIODIC */
            break;
 	  case shearing:
 	    *index = *(index - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
 	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (y-right).\n",
		    btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType!\n");
          }

        }

#else
 
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
 	    *index = *(index - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0]);
 	    break;
	  case BoundaryUndefined:
            break;
	  default:
	    ENZO_VFAIL("BoundaryType %"ISYM" not recognized (y-right).\n",
		    BoundaryType[field][1][1][bindex])
	  }

	}

#endif

  }
 
  /* set z inner (left) face */
 
  Sign = 1;
  if (FieldType == Velocity3) Sign = -1;
 
  if (BoundaryDimension[2] > 1 && GridOffset[2] == 0) {

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 2;
    face = 0;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (k = 0; k < StartIndex[2]; k++)
      for (i = 0; i < GridDims[0]; i++)
        for (j = 0; j < GridDims[1]; j++) {

          index = Field + i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
          bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "zl %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index + (2*StartIndex[2]-1 - 2*k)*GridDims[0]*GridDims[1]));
            break;
          case outflow:
            *index =       *(index + (  StartIndex[2]   - k)*GridDims[0]*GridDims[1]) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
            break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
#endif /* USE_PERIODIC */
            break;
 	  case shearing:
 	    *index = *(index + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
 	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (z-left).\n",btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType");
          }
        }

#else
 
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
 	    *index = *(index + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
 	    break;
	  case BoundaryUndefined:
            break;
	  default:
	    ENZO_VFAIL("BoundaryType %"ISYM" not recognized (z-left).\n",
		    BoundaryType[field][2][0][bindex])
	  }

	}

#endif

  }
 
  if (BoundaryDimension[2] > 1 && GridOffset[2]+GridDims[2] == BoundaryDimension[2]) {
 
    /* set z outer (right) face */

#ifdef OOC_BOUNDARY
    field = ExternalBoundaryField;
    dim = 2;
    face = 1;
    slabsize = cubesize/BoundaryDimension[dim];

    if (ExternalBoundaryTypeIO) {
      if ( ! SimpleConstantBoundary ) {
        READ_BT(bt_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      }
    }

    if (ExternalBoundaryValueIO)
      READ_BV(bv_buffer, field, dim, face, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);

    for (k = 0; k < GridDims[2]-EndIndex[2]-1; k++)
      for (i = 0; i < GridDims[0]; i++)
        for (j = 0; j < GridDims[1]; j++) {

          index = Field + i + j*GridDims[0] +
            (k + EndIndex[2]+1)*GridDims[1]*GridDims[0];
          bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];

//          if (bt_buffer[bindex] != 3) fprintf(stderr, "zr %"ISYM" %"ISYM"\n", bindex, bt_buffer[bindex]);
//          bt_buffer[bindex] = 3;

          if (SimpleConstantBoundary) {
            btb = periodic;
          } else {
            btb = bt_buffer[bindex];
          }

          switch (btb) {
          case reflecting:
            *index = Sign*(*(index - (2*k + 1)*GridDims[0]*GridDims[1]));
            break;
          case outflow:
            *index =       *(index + (-1 - k)*GridDims[0]*GridDims[1]) ;
            break;
          case inflow:
            *index = bv_buffer[bindex];
            break;
          case periodic:
#ifdef USE_PERIODIC
            *index = *(index - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
#endif /* USE_PERIODIC */
            break;
 	  case shearing:
 	    *index = *(index - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
 	    break;
          case BoundaryUndefined:
            break;
          default:
            fprintf(stderr, "IO BoundaryType %"ISYM" not recognized (z-right).\n",btb);
            fprintf(stderr, "field %"ISYM" dim %"ISYM" face %"ISYM" slab %"ISYM" bindex %"ISYM" btb %"ISYM"\n",
              field, dim, face, slabsize, bindex, bt_buffer[bindex]);
            ENZO_FAIL("Unrecognized IO BoundaryType");
          }

        }

#else
 
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
 	    *index = *(index - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1]);
 	    break;
 	  case BoundaryUndefined:
            break;
	  default:
	    fprintf(stderr, "BoundaryType %"ISYM" not recognized (z-right).\n",
		    BoundaryType[field][2][1][bindex]);
            fprintf(stderr, "field %"ISYM" dim %"ISYM"\n",field, dim);

	    ENZO_FAIL("Unrecognized IO BoundaryType!\n");
	  }

	}

#endif

  }

#ifdef OOC_BOUNDARY
  if (ExternalBoundaryTypeIO) {
    if ( ! SimpleConstantBoundary ) {
      delete [] bt_buffer;
    }
  }

  if (ExternalBoundaryValueIO)

    delete [] bv_buffer;
#endif
 
  return SUCCESS;
 
}
