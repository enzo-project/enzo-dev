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
#ifdef ECUDA
#include "cuPPM.h"
#endif


int grid::SolvePPM_DE(int CycleNumber, int NumberOfSubgrids, 
		      fluxes *SubgridFluxes[], float *CellWidthTemp[], 
		      Elong_int GridGlobalStart[], int GravityOn, 
		      int NumberOfColours, int colnum[])
{

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  



  int nxz, nyz, nzz, ixyz;
  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  ixyz = CycleNumber % GridRank;

  // compute pressure
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  float *Pressure = new float[size];
  this->ComputePressure(Time, Pressure);

#ifdef ECUDA
  cuPPMParameter PPMPara;
  cuPPMData PPMData;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);
  if (UseCUDA) {
    cuPPMInitParameter(&PPMPara, DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, GENum,
                       MAX_NUMBER_OF_BARYON_FIELDS,
                       GridGlobalStart, GridStartIndex, GridEndIndex,
                       GridDimension, CellWidthTemp[0][0], CellWidthTemp[1][0],
                       CellWidthTemp[2][0], Gamma, 
                       PPMDiffusionParameter, PPMFlatteningParameter,
                       PPMSteepeningParameter, ConservativeReconstruction,
                       PositiveReconstruction, GravityOn,
                       DualEnergyFormalism, DualEnergyFormalismEta1,
                       DualEnergyFormalismEta2, PressureFree, tiny_number,
                       RiemannSolver, RiemannSolverFallback, 
                       NumberOfColours, colnum);
    cuPPMInitData(&PPMData, PPMPara);
    cuPPMSetBaryon(PPMData, PPMPara, BaryonField, Pressure, AccelerationField);
  }
#endif

  int i,j,k,n;
  for (n = ixyz; n < ixyz+GridRank; n++) {

    // Update in x-direction
    if ((n % GridRank == 0) && nxz > 1) {
      if (UseCUDA == 0) 
	for (k = 0; k < GridDimension[2]; k++) {
	  if (this->xEulerSweep(k, NumberOfSubgrids, SubgridFluxes, 
				GridGlobalStart, CellWidthTemp, GravityOn, 
				NumberOfColours, colnum, Pressure) == FAIL) {
	    ENZO_VFAIL("Error in xEulerSweep.  k = %d\n", k)
	      }
	} // ENDFOR k
      else {
#ifdef ECUDA
        cuPPMSweep(PPMData, PPMPara, dtFixed, 0);
        cuPPMSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 
                               GridGlobalStart, PPMData, PPMPara, 0);
#endif
      }
    } // ENDIF x-direction

    // Update in y-direction
    if ((n % GridRank == 1) && nyz > 1) {
      if (UseCUDA == 0) 
	for (i = 0; i < GridDimension[0]; i++) {
	  if (this->yEulerSweep(i, NumberOfSubgrids, SubgridFluxes, 
				GridGlobalStart, CellWidthTemp, GravityOn, 
				NumberOfColours, colnum, Pressure) == FAIL) {
	    ENZO_VFAIL("Error in yEulerSweep.  i = %d\n", i)
	      }
	} // ENDFOR i
      else {
#ifdef ECUDA
        cuPPMSweep(PPMData, PPMPara, dtFixed, 1);
        cuPPMSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 
                               GridGlobalStart, PPMData, PPMPara, 1);
#endif
      }
    } // ENDIF y-direction
      
      // Update in z-direction
    if ((n % GridRank == 2) && nzz > 1) {
      if (UseCUDA == 0) 
	for (j = 0; j < GridDimension[1]; j++) {
	  if (this->zEulerSweep(j, NumberOfSubgrids, SubgridFluxes, 
				GridGlobalStart, CellWidthTemp, GravityOn, 
				NumberOfColours, colnum, Pressure) == FAIL) {
	    ENZO_VFAIL("Error in zEulerSweep.  j = %d\n", j)

	      }
	} // ENDFOR j
      else {
#ifdef ECUDA
	cuPPMSweep(PPMData, PPMPara, dtFixed, 2);
	cuPPMSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 
			       GridGlobalStart, PPMData, PPMPara, 2);
#endif
      }
    } // ENDIF z-direction
      
  } // ENDFOR n

#ifdef ECUDA
  if (UseCUDA != 0) {
    cuPPMGetBaryon(PPMData, PPMPara, BaryonField);
    cuPPMDestroy(PPMData, PPMPara);
  }
#endif

  // make energies consistent for polytropic equations of state
  if (EOSType > 0) 
    this->ComputePressure(Time, Pressure);


  delete [] Pressure;

  return SUCCESS;

}
