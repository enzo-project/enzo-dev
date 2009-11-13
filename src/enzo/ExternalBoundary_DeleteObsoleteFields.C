/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (DELETES OBSOLETE EXTERNAL BOUNDARY VALUES)
/
/  written by: John Wise
/  date:       November, 2009
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

int ExternalBoundary::DeleteObsoleteFields(int *ObsoleteFields,
					   int NumberOfObsoleteFields)
{

  int i, j, k, dim, field;
  int NumberOfDeletions = 0;

  for (field = 0; field < NumberOfObsoleteFields; field++)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (BoundaryFieldType[i] == ObsoleteFields[field]) {

	/* Delete values and shift BoundaryFieldType and BoundaryType
	   back */

	for (j = i; j < MAX_NUMBER_OF_BARYON_FIELDS-1; j++)
	  BoundaryFieldType[j] = BoundaryFieldType[j+1];
	BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS-1] = FieldUndefined;

	for (dim = 0; dim < BoundaryRank; dim++)
	  if (BoundaryDimension[dim] > 1) {
	    for (j = 0; j < 2; j++) {

	      delete [] BoundaryType[i][dim][j];
	      for (k = i; k < MAX_NUMBER_OF_BARYON_FIELDS-1; k++)
		BoundaryType[k][dim][j] = BoundaryType[k+1][dim][j];
	      BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS-1][dim][j] = NULL;

	    } // ENDFOR direction (face)
	  } // ENDIF dimension > 1

	NumberOfDeletions++;

	break;

      } // ENDIF matching field

  NumberOfBaryonFields -= NumberOfDeletions;

  return SUCCESS;

}
