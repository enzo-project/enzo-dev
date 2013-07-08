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
#include "ErrorExceptions.h"
#include "performance.h"
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
  if (!(RandomForcing)) return SUCCESS;
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  int i, dim;

  LCAPERF_START("grid_AddRandomForcing");

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("GARF: Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* error check. */
 
  if (RandomForcingField[0] == NULL)
    ERROR_MESSAGE;
 
  int corneri=  GridStartIndex[0] + GridDimension[0]*(GridStartIndex[1]+GridStartIndex[2]*GridDimension[1]);
  if (RandomForcingField[0][0] == 0.0)
    ERROR_MESSAGE;
  fprintf(stderr, "TopGridTimeStep: %g\n", dtTopGrid);
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
 
  if ( HydroMethod == MHD_Li )
      MHDCT_ConvertEnergyToSpecificC();//See docs or Grid_MHDCTEnergyToggle.C for if/when this is done

  /* do not do the update if using ZEUS */
 
  if (HydroMethod != Zeus_Hydro && EquationOfState == 0)
    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++){
	BaryonField[TENum][i] +=
	  BaryonField[Vel1Num+dim][i]*RandomForcingField[dim][i]*levelNorm +
	  0.5*RandomForcingField[dim][i]*levelNorm*RandomForcingField[dim][i]*levelNorm;
      }
 
  /* add velocity perturbation to velocity fields. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
	BaryonField[Vel1Num+dim][i] += RandomForcingField[dim][i]*levelNorm;

  if ( HydroMethod == MHD_Li )
      MHDCT_ConvertEnergyToConservedC();  //See docs or Grid_MHDCTEnergyToggle.C for if/when this is done
 
  LCAPERF_STOP("grid_AddRandomForcing");
  return SUCCESS;
 
}
