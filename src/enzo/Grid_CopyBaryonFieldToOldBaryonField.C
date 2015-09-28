/***********************************************************************
/
/  GRID CLASS (COPY BARYON FIELD TO OLD BARYON FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness / Brian O'Shea
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
#include "performance.h"
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
 
  LCAPERF_START("grid_CopyBaryonFieldToOldBaryonField");

  /* compute the field size */
 
  int size = 1;

  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* copy fields */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
 
    /* Check to make sure BaryonField exists. */
 
    if (BaryonField[field] == NULL) {
      ENZO_FAIL("BaryonField missing.\n");
    }

    /* Create OldBaryonField if necessary. */
 
    if ((OldBaryonField[field] == NULL))
      OldBaryonField[field] = new float[size];
 
    /* Copy. */
 
    for (i = 0; i < size; i++)
      OldBaryonField[field][i] = BaryonField[field][i];
 
  } // end loop over fields

  if(UseMHDCT){   
    for(field=0;field<3;field++){

      if(MagneticField[field] == NULL )
	ENZO_FAIL("MagneticField mising in CopyBaryonFieldToOldBaryonField");
      
      if(OldMagneticField[field] == NULL) {
	OldMagneticField[field] = new float[MagneticSize[field]];
      }
      
      for(i=0;i<MagneticSize[field];i++){
	OldMagneticField[field][i] = MagneticField[field][i];
      }

      if(CenteredB[field] == NULL)
	ENZO_FAIL("CenteredB missing in CopyBaryonFieldToOldBaryonField");
      
      if(OldCenteredB[field] == NULL) {
	OldCenteredB[field] = new float[size];
      }
      
      for(i=0;i<size;i++){
	OldCenteredB[field][i] = CenteredB[field][i];
      }

    }//for(field < 3;)
  }//end if(UseMHDCT)

  // AccelerationHack

#ifdef SAB

  // Mod from Brian O'Shea, 8th August 2006
  // In case there are no baryon fields

  if( (SelfGravity || UniformGravity || PointSourceGravity || DiskGravity ) && (NumberOfBaryonFields > 0) ) {

    for(field = 0; field < GridRank; field++) {

      if(AccelerationField[field] != NULL) {
        if( OldAccelerationField[field] == NULL ) {
          OldAccelerationField[field] = new float[size];
        }
        for(i=0;i<size;i++) {
          OldAccelerationField[field][i] = AccelerationField[field][i];
        }

      }else{

        ENZO_FAIL("Error-- in CopyBF to Old, no AccelerationField.\n");


      }

    }  //field, GravityFlags.

  }

#endif /* SAB */
 
  LCAPERF_STOP("grid_CopyBaryonFieldToOldBaryonField");
  return SUCCESS;
 
}
