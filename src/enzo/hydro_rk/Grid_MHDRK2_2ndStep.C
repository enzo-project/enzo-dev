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
int MHDTimeUpdate_CUDA(float **Prim, int GridDimension[], 
			int GridStartIndex[], int GridEndIndex[], int GridRank,
		        float dtdx, float dt, float C_h, float C_p);

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


  double time1 = ReturnWallTime();

  float *Prim[NEQ_MHD+NSpecies+NColor];

  this->ReturnHydroRKPointers(Prim,0);


#ifdef ECUDA
  if (UseCUDA == 1) {
    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time3 = ReturnWallTime();
    if (MHDTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
			    dtdx, dtFixed, C_h, C_p) == FAIL) {
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

    for (int field = 0; field < NEQ_MHD; field++) {
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

  } // if (UseCUDA)
#endif // ifdef ECUDA

  if (StellarWindFeedback)
    this->ReduceWindBoundary();

  /* Compute dU */

  float *dU[NEQ_MHD+NSpecies+NColor];
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
  }

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

  return SUCCESS;

}
