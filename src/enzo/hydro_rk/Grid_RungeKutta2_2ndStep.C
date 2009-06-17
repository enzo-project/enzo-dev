/***********************************************************************
/
/  GRID CLASS (SECOND STEP OF RUNGE-KUTTA INTEGRATION)
/
/  written by: Peng Wang
/  date:       May, 2007
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


int grid::RungeKutta2_2ndStep(int CycleNumber, fluxes *SubgridFluxes[], 
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
    B1Num, B2Num, B3Num, HMNum, H2INum, H2IINum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  double time1 = ReturnWallTime();

  float *dU[NEQ_HYDRO+NSpecies+NColor];
  float *Prim[NEQ_HYDRO+NSpecies+NColor];

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);


  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
    for (int i = 0; i < activesize; i++) {
      dU[field][i] = 0.0;
    }
  }

  Prim[iden ] = BaryonField[DensNum];
  Prim[ivx  ] = BaryonField[Vel1Num];
  Prim[ivy  ] = BaryonField[Vel2Num];
  Prim[ivz  ] = BaryonField[Vel3Num];
  Prim[ietot] = BaryonField[TENum];
  if (DualEnergyFormalism) {
    Prim[ieint] = BaryonField[GENum];
  }


  // copy species field
  for (int ns = NEQ_HYDRO; ns < NEQ_HYDRO+NSpecies; ns++) {
    // change species from density to mass fraction
    for (int n = 0; n < size; n++) {
      BaryonField[ns][n] /= BaryonField[iden][n];
    }
    Prim[ns] = BaryonField[ns];
  }

  // copy color field
  for (int nc = NEQ_HYDRO+NSpecies; nc < NEQ_HYDRO+NSpecies+NColor; nc++) {
    Prim[nc] = BaryonField[nc];
  }

  int fallback = 0;

  if (this->Hydro3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 
		    0.5, fallback) == FAIL) {
    return FAIL;
  }

  this->SourceTerms(dU);

  if (this->UpdatePrim(dU, 0.5, 0.5) == FAIL) {
    // fall back to zero order scheme
    printf("Falling back to zero order at RK 2nd step\n");

    this->CopyOldBaryonFieldToBaryonField();
    // change species from density to mass fraction
    for (int ns = NEQ_HYDRO; ns < NEQ_HYDRO+NSpecies; ns++) {
      for (int n = 0; n < size; n++) {
	BaryonField[ns][n] /= BaryonField[iden][n];
      }
    }

    this->ZeroFluxes(SubgridFluxes, NumberOfSubgrids);
    fallback = 1;
    if (this->Hydro3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids,
		      0.5, fallback) == FAIL) {
      return FAIL;
    }
    this->SourceTerms(dU);
    if (this->UpdatePrim(dU, 0.5, 0.5) == FAIL) {
      printf("Fallback failed, give up...\n");
      return FAIL;
    }
    return FAIL;
  }

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    delete [] dU[field];
  }

  //  PerformanceTimers[1] += ReturnWallTime() - time1;

  return SUCCESS;

}
