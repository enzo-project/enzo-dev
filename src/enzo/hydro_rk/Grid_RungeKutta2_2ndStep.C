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
int HydroTimeUpdate_CUDA(float **Prim, int GridDimension[], int GridStartIndex[], int GridEndIndex[], int GridRank,
		          float dtdx, float dt);


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
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  double time1 = ReturnWallTime();

  float *Prim[NEQ_HYDRO+NSpecies+NColor];

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);



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

#ifdef ECUDA
  if (UseCUDA == 1) {

    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time3 = ReturnWallTime();
    if (HydroTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
		      	      dtdx, dtFixed) == FAIL) {
      printf("RK2: MHDTimeUpdate_CUDA failed.\n");
      return FAIL;
    }
    
    double time2 = ReturnWallTime();

    for (int field = ivx; field <= ietot; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid =i + (j + k*GridDimension[1])*GridDimension[0];
	    BaryonField[field][igrid] *= BaryonField[iden][igrid];
	    OldBaryonField[field][igrid] *= OldBaryonField[iden][igrid];
	  }
	}
      }
    }

    for (int field = 0; field < NEQ_HYDRO; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid =i + (j + k*GridDimension[1])*GridDimension[0];
	    BaryonField[field][igrid] = 0.5*(OldBaryonField[field][igrid] + BaryonField[field][igrid]);
	  }
	}
      }
    }

    for (int field = ivx; field <= ietot; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	    BaryonField[field][igrid] /= BaryonField[iden][igrid];
	    OldBaryonField[field][igrid] /= OldBaryonField[iden][igrid];
	  }
	}
      }
    }

    return SUCCESS;

  } // if UseCUDA == 1
#endif /* ECUDA */

  float *dU[NEQ_HYDRO+NSpecies+NColor];
  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
    for (int i = 0; i < activesize; i++) {
      dU[field][i] = 0.0;
    }
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
