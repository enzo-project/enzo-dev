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
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  int i, dim;
  float *ForcingField[GridRank];
  float levelNorm = (*norm)*dtFixed/dtTopGrid;

  // use original RandomForcing field
  if (RandomForcing) { 

    for (dim = 0; dim < GridRank; dim++)
        ForcingField[dim] = RandomForcingField[dim];
  
    if (dtTopGrid == 0.0)
      ERROR_MESSAGE;
    levelNorm = (*norm)*dtFixed/dtTopGrid;

  // use StochasticForcing in combination with MHDCT
  // for RK2 based solvers, the forcing happens in the Grid_*SourceTerms.C
  } else if (UseDrivingField && HydroMethod == MHD_Li) {
    int Drive1Num, Drive2Num, Drive3Num;
    if (this->IdentifyDrivingFields(Drive1Num, Drive2Num, Drive3Num) == FAIL) {
      printf("grid::AddRandomForcing: canot identify driving fields.\n");
      return FAIL;
    }

    for (dim = 0; dim < GridRank; dim++)
      ForcingField[dim] = BaryonField[Drive1Num+dim];

    // this ensures, that a) norm is set and b) we'll always use the proper 
    // timestep. This may change when StochasticForcing has been tested for
    // other configurations than static, uniform grids
    levelNorm = dtFixed;

  /* Return if this doesn't concern us. */
  } else  
    return SUCCESS;


  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;


  LCAPERF_START("grid_AddRandomForcing");

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("GARF: Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* error check. */
 
  if (ForcingField[0] == NULL)
    ERROR_MESSAGE;
 
  if (ForcingField[0][0] == 0.0)
    ERROR_MESSAGE;
  //fprintf(stderr, "TopGridTimeStep: %g\n", dtTopGrid);
 
  /* compute the field size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* update total energy first. */
 
 
  if (levelNorm <= 0.0)
    WARNING_MESSAGE;
 
  if ( HydroMethod == MHD_Li )
      MHDCT_ConvertEnergyToSpecificC();//See docs or Grid_MHDCTEnergyToggle.C for if/when this is done

  /* do not do the update if using ZEUS */
 
  if (HydroMethod != Zeus_Hydro && EquationOfState == 0)
    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++){
	BaryonField[TENum][i] +=
	  BaryonField[Vel1Num+dim][i]*ForcingField[dim][i]*levelNorm +
	  0.5*ForcingField[dim][i]*levelNorm*ForcingField[dim][i]*levelNorm;
      }
 
  /* add velocity perturbation to velocity fields. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
	BaryonField[Vel1Num+dim][i] += ForcingField[dim][i]*levelNorm;

  if ( HydroMethod == MHD_Li )
      MHDCT_ConvertEnergyToConservedC();  //See docs or Grid_MHDCTEnergyToggle.C for if/when this is done
 
  LCAPERF_STOP("grid_AddRandomForcing");
  return SUCCESS;
 
}
