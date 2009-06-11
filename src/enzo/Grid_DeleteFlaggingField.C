/***********************************************************************
/
/  GRID CLASS (DELETE THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
// Delete the flagging field
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::DeleteFlaggingField()
{
  if (FlaggingField == NULL)
    fprintf(stderr, "DeleteFlaggingField: Warning, field not present.\n");
  else
    delete [] FlaggingField;
  FlaggingField = NULL;
}
