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
    if ((n % 3 == 0) && nxz > 1) {
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	if (this->xEulerSweep(k, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  fprintf(stderr, "Error in xEulerSweep.  k = %d\n", k);
	  ENZO_FAIL("");
	}
      } // ENDFOR k
    } // ENDIF x-direction

    // Update in y-direction
    if ((n % 3 == 1) && nyz > 1) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	if (this->yEulerSweep(i, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  fprintf(stderr, "Error in yEulerSweep.  i = %d\n", i);
	  ENZO_FAIL("");
	}
      } // ENDFOR i
    } // ENDIF y-direction

    // Update in z-direction
    if ((n % 3 == 2) && nzz > 1) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	if (this->zEulerSweep(j, NumberOfSubgrids, SubgridFluxes, 
			      GridGlobalStart, CellWidthTemp, GravityOn, 
			      NumberOfColours, colnum) == FAIL) {
	  fprintf(stderr, "Error in zEulerSweep.  j = %d\n", j);
	  ENZO_FAIL("");
	}
      } // ENDFOR j
    } // ENDIF z-direction

  } // ENDFOR n

  return SUCCESS;

}
