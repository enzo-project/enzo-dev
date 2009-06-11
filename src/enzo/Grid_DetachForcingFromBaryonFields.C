/***********************************************************************
/
/  GRID CLASS (REMOVE RANDOM FORCING FIELDS FROM THE LIST OF BARYON FIELDS)
/
/  written by: Alexei Kritsuk
/  date:       January 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
 
int grid::DetachForcingFromBaryonFields()
{
 
  /* Detach RandomForcingFields from BaryonFields. */
 
  for (int dim = 0; dim < GridRank; dim++) {
    BaryonField[NumberOfBaryonFields-dim-1] = NULL;
    FieldType[NumberOfBaryonFields-dim-1] = FieldUndefined;
   }
  NumberOfBaryonFields -= GridRank;
 
  //  if (debug) printf("ForcingDetached[%"ISYM"]\n", ProcessorNumber);
 
  return SUCCESS;
 
}
