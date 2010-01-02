/***********************************************************************
/
/  GRID CLASS (DELETE OBSOLETE FIELDS WHEN RESTARTING)
/
/  written by: John Wise
/  date:       November, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

bool FirstTime = true;

int grid::DeleteObsoleteFields(int *ObsoleteFields, 
			       int NumberOfObsoleteFields)
{

  int i, j, field, NumberOfDeletions = 0;

  for (field = 0; field < NumberOfObsoleteFields; field++)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (FieldType[i] == ObsoleteFields[field]) {

	if (debug && FirstTime)
	  printf("Deleting unused field %d (FieldType = %d = %s)\n", 
		 i, ObsoleteFields[field], DataLabel[i]);

	/* Delete field */

	if (MyProcessorNumber == ProcessorNumber)
	  delete [] BaryonField[i];

	/* Shift FieldType and BaryonField back */

	for (j = i; j < MAX_NUMBER_OF_BARYON_FIELDS-1; j++)
	  FieldType[j] = FieldType[j+1];
	FieldType[MAX_NUMBER_OF_BARYON_FIELDS-1] = FieldUndefined;

	if (MyProcessorNumber == ProcessorNumber) {
	  for (j = i; j < MAX_NUMBER_OF_BARYON_FIELDS-1; j++)
	    BaryonField[j] = BaryonField[j+1];
	  BaryonField[MAX_NUMBER_OF_BARYON_FIELDS-1] = NULL;
	}

	if (FirstTime) {
	  for (j = i; j < MAX_NUMBER_OF_BARYON_FIELDS-1; j++) {
	    DataLabel[j] = DataLabel[j+1];
	    DataUnits[j] = DataUnits[j+1];
	  }
	  DataLabel[MAX_NUMBER_OF_BARYON_FIELDS-1] = NULL;
	  DataUnits[MAX_NUMBER_OF_BARYON_FIELDS-1] = NULL;
	}

	NumberOfDeletions++;
	break;

      } // ENDIF matching field

  NumberOfBaryonFields -= NumberOfDeletions;
  FirstTime = false;

  return SUCCESS;

}
