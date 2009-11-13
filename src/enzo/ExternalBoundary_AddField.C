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

int ExternalBoundary::AddField(int FieldType)
{

  int i, j, dim, ifield, size;

  ifield = NumberOfBaryonFields;
  BoundaryFieldType[ifield] = FieldType;

  /* loop over faces */

  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] > 1) {

      /* Calculate boundary size */

      size = 1;
      for (i = 0; i < BoundaryRank; i++)
	if (i != dim)
	  size *= BoundaryDimension[i];

      for (i = 0; i < 2; i++) {

	/* allocate room for BoundaryType */

	BoundaryType[ifield][dim][i] = new boundary_type[size];

	/* assign boundary of the new field the same as the first field */

	for (j = 0; j < size; j++)
	  BoundaryType[ifield][dim][i][j] = BoundaryType[0][dim][i][j];

      }   // end of loop over dims

    } // ENDIF BoundaryDimension > 1

  NumberOfBaryonFields++;

  return SUCCESS;

}
