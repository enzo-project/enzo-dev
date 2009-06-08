/***********************************************************************
/
/  GRID CLASS (ADD RANDOM FORCING FIELDS TO VELOCITIES)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::AddRandomForcing(float * norm, float dtTopGrid)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  int i, dim;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "GARF: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
 
  /* error check. */
 
  if (RandomForcingField[0] == NULL)
    ERROR_MESSAGE;
 
  if (RandomForcingField[0][0] == 0.0)
    ERROR_MESSAGE;
 
  if (dtTopGrid == 0.0)
    ERROR_MESSAGE;
 
 
  /* compute the field size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* update total energy first. */
 
 
  float levelNorm = (*norm)*dtFixed/dtTopGrid;
  if (levelNorm <= 0.0)
    WARNING_MESSAGE;
 
  /* do not do the update if using ZEUS */
 
  if (HydroMethod != Zeus_Hydro)
    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++)
	BaryonField[TENum][i] +=
	  BaryonField[Vel1Num+dim][i]*
	  RandomForcingField[dim][i]*levelNorm +
	  0.5*RandomForcingField[dim][i]*levelNorm*
	  RandomForcingField[dim][i]*levelNorm;
 
  /* add velocity perturbation to velocity fields. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[Vel1Num+dim][i] += RandomForcingField[dim][i]*levelNorm;
 
  return SUCCESS;
 
}
