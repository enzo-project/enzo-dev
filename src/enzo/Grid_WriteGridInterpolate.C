/***********************************************************************
/
/  GRID CLASS (INTERPOLATE AND OUTPUT GRID TO ARBITRARY TIME)
/
/  written by: Greg Bryan
/  date:       April, 2000
/  modified1:
/
/  PURPOSE:  This routine interpolates grid data to the time passed
/            in and then call the regular grid i/o routine.  It is
/            intended for outputs at arbitrary times, rather than just
/            at the end of a step, as the regular write grid routine
/            assumes.
/
/  RETURNS:
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
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
/* This macro converts a float and writes it to the local buffer, which,
   when full is written to the file pointed to by fptr. */
 
 
int grid::WriteGridInterpolate(FLOAT WriteTime, FILE *fptr, char *base_name,
			       int grid_id)
{
 
  int dim, i, field, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute coefficient factors for linear interpolation in time.
     Note: interp = coef1*Old + coef2*New. */
 
  float coef1 = 0, coef2 = 1;
  if (Time != WriteTime) {
    if (Time <= OldTime) {
      ENZO_FAIL("WGI: fields are at the same time or worse.\n");
    } else {
      coef1 = max((Time - WriteTime)/
		  (Time - OldTime), 0.0);
      coef2  = (1.0 - coef1);
    }
  }
 
  /* Interpolate grid to given time and save old variables. */
 
  float *SavedBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  if (coef2 != 1 && MyProcessorNumber == ProcessorNumber)
    for (field = 0; field < NumberOfBaryonFields; field++) {
      SavedBaryonField[field] = BaryonField[field];
      BaryonField[field] = new float[size];
      for (i = 0; i < size; i++)
	BaryonField[field][i] = coef1*OldBaryonField[field][i] +
	                        coef2*SavedBaryonField[field][i];
    }
 
  /* Move particles to given time. */
 
  float TimeDifference = WriteTime - Time;
  if (this->UpdateParticlePosition(TimeDifference) == FAIL) {
    ENZO_FAIL("Error in grid->UpdateParticlePosition.\n");
  }
 
  /* Write grid (temporarily replace Time with WriteTime). */
 
  FLOAT SavedTime = Time;
  Time = WriteTime;
  if (this->WriteGrid(fptr, base_name, grid_id) == FAIL) {
    ENZO_FAIL("Error in grid->WriteGrid.\n");
  }
  Time = SavedTime;
 
  /* Move particles back. */
 
  this->UpdateParticlePosition(-TimeDifference);
 
  /* Copy saved baryon fields back. */
 
  if (coef2 != 1 && MyProcessorNumber == ProcessorNumber)

    for (field = 0; field < NumberOfBaryonFields; field++) {
      delete [] BaryonField[field];
      BaryonField[field] = SavedBaryonField[field];
    }
 
  return SUCCESS;
}
 
