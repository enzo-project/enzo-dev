/***********************************************************************
/
/  GRID CLASS (SUM MASS FLAGGING FIELD FROM ALL PROCESSORS)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  
/
/  PURPOSE:
/     When we're finding new subgrids, we keep the particles from subgrids
/       on their original processor to distribute memory usage.  Because
/       of this, we have to sum the mass flagging field before we look
/       for flagged cells.
/
/  NOTE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif 
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
 
int grid::CollectParticleMassFlaggingField(void)
{
 
  /* If serial, do nothing. */
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
  /* error check */
 
  if (ParticleMassFlaggingField == NULL) {
    ENZO_FAIL("ParticleMassFlaggingField is undefined.\n");
  }

  /* compute size */

  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Collect the sum on the grid's main processor */

#ifdef USE_MPI

  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count = size;

  float *buffer = new float[size];
  for (i = 0; i < size; i++)
    buffer[i] = ParticleMassFlaggingField[i];

  MPI_Reduce(buffer, ParticleMassFlaggingField, Count, DataType, MPI_SUM, 
	     ProcessorNumber, MPI_COMM_WORLD);

  delete [] buffer;

#endif /* USE_MPI */  

  /* Delete MassFlaggingField if this isn't the main processor */

  if (MyProcessorNumber != ProcessorNumber) {

    delete [] ParticleMassFlaggingField;
    ParticleMassFlaggingField = NULL;
  }
 
  return SUCCESS;
}
