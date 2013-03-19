/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR MHD Li (MHDCT) method
/
/
/  written by: dcollins
/  date:       March 19, 2013, 3:54 pm.  
/  modified1:
/
/  PURPOSE:  
/
/  RETURNS:
/    SUCCESS or FAIL
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

int grid::SolveMHD_Li(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[])
{

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  ENZO_FAIL("MHDLi called.");
  
  return SUCCESS;
}
