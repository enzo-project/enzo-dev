/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR EULER SOLVER)
/
/  written by: John Wise
/  date:       May, 2007
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//


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
//#include "fortran.def"

int grid::SolvePPM_DE(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[])
{

  int nxz, nyz, nzz, ixyz;

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  ixyz = CycleNumber % GridRank;

  int i,j,k,n;
  for (n = ixyz; n < ixyz+GridRank; n++) {

    // Update in x-direction
    if ((n % GridRank == 0) && nxz > 1) {
      for (k = 0; k < GridDimension[2]; k++) {
	if (this->xEulerSweep(k, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  ENZO_VFAIL("Error in xEulerSweep.  k = %d\n", k)
	}
      } // ENDFOR k
    } // ENDIF x-direction

    // Update in y-direction
    if ((n % GridRank == 1) && nyz > 1) {
      for (i = 0; i < GridDimension[0]; i++) {
	if (this->yEulerSweep(i, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  ENZO_VFAIL("Error in yEulerSweep.  i = %d\n", i)
	}
      } // ENDFOR i
    } // ENDIF y-direction

    // Update in z-direction
    if ((n % GridRank == 2) && nzz > 1) {
      for (j = 0; j < GridDimension[1]; j++) {
	if (this->zEulerSweep(j, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  ENZO_VFAIL("Error in zEulerSweep.  j = %d\n", j)

	}
      } // ENDFOR j
    } // ENDIF z-direction

  } // ENDFOR n

  return SUCCESS;

}
