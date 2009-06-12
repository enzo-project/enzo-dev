/***********************************************************************
/
/  GRID CLASS (RUNGE-KUTTA SECOND STEP)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"


double ReturnWallTime();


int grid::MHDRK2_2ndStep(int CycleNumber, fluxes *SubgridFluxes[], 
			      int NumberOfSubgrids, int level,
			      ExternalBoundary *Exterior)
  /*
    NumberOfSubgrids: the actual number of subgrids + 1
    SubgridFluxes[NumberOfSubgrids]
  */
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (NumberOfBaryonFields == 0) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, PhiNum, HMNum, H2INum, H2IINum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  double time1 = ReturnWallTime();

  float *dU[NEQ_MHD+NSpecies+NColor];
  float *Prim[NEQ_MHD+NSpecies+NColor];

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }
  
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);
  }


  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
  }

  if (StellarWindFeedback) {
    this->ReduceWindBoundary();
  }

  Prim[iden ] = BaryonField[DensNum];
  Prim[ivx  ] = BaryonField[Vel1Num];
  Prim[ivy  ] = BaryonField[Vel2Num];
  Prim[ivz  ] = BaryonField[Vel3Num];
  Prim[ietot] = BaryonField[TENum];
  if (DualEnergyFormalism) {
    Prim[ieint] = BaryonField[GENum];
  }

  Prim[iBx  ] = BaryonField[B1Num];
  Prim[iBy  ] = BaryonField[B2Num];
  Prim[iBz  ] = BaryonField[B3Num];
  Prim[iPhi ] = BaryonField[PhiNum];

  /* Copy species field */

  for (int ns = NEQ_MHD; ns < NEQ_MHD+NSpecies; ns++) {
    /* change species from density to mass fraction */
    for (int n = 0; n < size; n++) {
      BaryonField[ns][n] /= BaryonField[iden][n];
    }
    Prim[ns] = BaryonField[ns];
  }

  /* Copy color field */

  for (int nc = NEQ_MHD+NSpecies; nc < NEQ_MHD+NSpecies+NColor; nc++) {
    Prim[nc] = BaryonField[nc];
  }

  /* Compute dU */

  int fallback = 0;
  if (this->MHD3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 
		  0.5, fallback) == FAIL) {
    return FAIL;
  }

  /* Add source terms */

  this->MHDSourceTerms(dU);

  /* Update primitive variables */

  if (this->UpdateMHDPrim(dU, 0.5, 0.5) == FAIL) {
    return FAIL;
  }

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] dU[field];
  }

  //  PerformanceTimers[1] += ReturnWallTime() - time1;

  return SUCCESS;

}
