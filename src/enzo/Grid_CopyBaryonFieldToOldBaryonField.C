/***********************************************************************
/
/  GRID CLASS (COPY BARYON FIELD TO OLD BARYON FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness / Biran O'Shea
/  date:       4th June 2006
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Copy the current baryon fields to the old baryon fields
//   (allocate old baryon fields if they don't exist).
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::CopyBaryonFieldToOldBaryonField()
{

  int i, field;
 
  /* update the old baryon field time */
 
  OldTime = Time;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* compute the field size */
 
  int size = 1;

/* BUG reported 4th June 2006
  for (int dim = 0; dim < GridDimension[dim]; dim++)
    size *= GridDimension[dim];
*/

  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* copy fields */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
 
    /* Check to make sure BaryonField exists. */
 
    if (BaryonField[field] == NULL) {
      fprintf(stderr, "BaryonField missing.\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
 
    /* Create OldBaryonField if necessary. */
 
    if ((OldBaryonField[field] == NULL))
      OldBaryonField[field] = new float[size];
 
    /* Copy. */
 
    for (i = 0; i < size; i++)
      OldBaryonField[field][i] = BaryonField[field][i];
 
  } // end loop over fields

  // AccelerationHack

#ifdef SAB

  // Mod from Brian O'Shea, 8th August 2006
  // In case there are no baryon fields

  if( (SelfGravity || UniformGravity || PointSourceGravity) && (NumberOfBaryonFields > 0) ) {

    for(field = 0; field < GridRank; field++) {

      if(AccelerationField[field] != NULL) {
        if( OldAccelerationField[field] == NULL ) {
          OldAccelerationField[field] = new float[size];
        }
        for(i=0;i<size;i++) {
          OldAccelerationField[field][i] = AccelerationField[field][i];
        }

      }else{

        fprintf(stderr,"Error-- in CopyBF to Old, no AccelerationField.\n");
        ENZO_FAIL("Error in: "__FILE__);

      }

    }  //field, GravityFlags.

  }

#endif /* SAB */
 
  return SUCCESS;
 
}
