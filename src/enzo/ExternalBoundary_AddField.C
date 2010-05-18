/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (ADD EXTERNAL BOUNDARY VALUES FOR A NEW FIELD)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// This routine reads the external boundary from the provided file pointer
//

/* function prototypes */

int READ_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, 
	    int BoundaryDimension[], int BoundaryRank, int Nfields);
int WRITE_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, 
	     int BoundaryDimension[], int BoundaryRank, int NumberOfBaryonFields);

int ExternalBoundary::AddField(int FieldType)
{

  int i, j, dim, field, ifield, size;
#ifdef OOC_BOUNDARY
  int slabsize;
  boundary_type *bt_buffer0, *bt_buffer1;
#endif

  ifield = NumberOfBaryonFields;
  BoundaryFieldType[ifield] = FieldType;

#ifdef OOC_BOUNDARY
  if (SimpleConstantBoundary) {
    NumberOfBaryonFields++;
    return SUCCESS;
  }
#endif

  /* loop over faces */

  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] > 1) {

      /* Calculate boundary size */

      size = 1;
      for (i = 0; i < BoundaryRank; i++)
	if (i != dim)
	  size *= BoundaryDimension[i];

#ifdef OOC_BOUNDARY

      /* Read the boundary types for the first field and duplicate
	 them for all baryon fields, including the added field.  This
	 won't work for simulations with different boundary types for
	 different fields.  Does this ever happen?  */

      slabsize = size/BoundaryDimension[dim];
      bt_buffer0  = new boundary_type[slabsize];
      bt_buffer1  = new boundary_type[slabsize];
      READ_BT(bt_buffer0, 0, dim, 0, slabsize, BoundaryDimension, BoundaryRank,
	      NumberOfBaryonFields);
      READ_BT(bt_buffer1, 0, dim, 1, slabsize, BoundaryDimension, BoundaryRank,
	      NumberOfBaryonFields);

      // Overwrite the boundary type file with the added field
      for (field = 0; field < NumberOfBaryonFields+1; field) {
	WRITE_BT(bt_buffer0, field, dim, 0, slabsize, BoundaryDimension,
		 BoundaryRank, NumberOfBaryonFields+1);
	WRITE_BT(bt_buffer1, field, dim, 1, slabsize, BoundaryDimension,
		 BoundaryRank, NumberOfBaryonFields+1);
      }
      delete [] bt_buffer0;
      delete [] bt_buffer1;
#else
      for (i = 0; i < 2; i++) {

	/* allocate room for BoundaryType */
	
	BoundaryType[ifield][dim][i] = new boundary_type[size];

	/* assign boundary of the new field the same as the first field */

	for (j = 0; j < size; j++)
	  BoundaryType[ifield][dim][i][j] = BoundaryType[0][dim][i][j];

      }   // end of loop over dims
#endif /* OOC_BOUNDARY */

    } // ENDIF BoundaryDimension > 1

  NumberOfBaryonFields++;

  return SUCCESS;

}
