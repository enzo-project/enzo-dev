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
#include "ErrorExceptions.h"
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
		        float dtdx, float dt, float C_h, float C_p, float cTheta_Limiter);

int grid::MHDRK2_2ndStep(fluxes *SubgridFluxes[], 
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
  float *OldPrim[NEQ_MHD+NSpecies+NColor];

#ifdef ECUDADEBUG
  printf("in Grid_MHDRK_2ndStep.C.\n");
  for (int j=30; j < 33; j++) 
    for (int i=0; i < 9; i++) printf("BaryonField[%"ISYM"][%"ISYM"] = %"GSYM" \n", i, j, BaryonField[i][j]);
#endif

  this->ReturnHydroRKPointers(Prim, false);  
  this->ReturnOldHydroRKPointers(OldPrim, false);  

#ifdef ECUDADEBUG
  printf("in Grid_MHDRK_2ndStep.C.\n");
  for (int j=30; j < 33; j++) 
    for (int i=0; i < 9; i++) printf("Prim[%"ISYM"][%"ISYM"] = %"GSYM" \n", i, j, Prim[i][j]);
#endif


#ifdef ECUDA
  if (UseCUDA == 1) {
    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time3 = ReturnWallTime();
    if (MHDTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
			   dtdx, dtFixed, C_h, C_p, Theta_Limiter) == FAIL) {
      printf("RK2: MHDTimeUpdate_CUDA failed.\n");
      return FAIL;
    }
    
    double time2 = ReturnWallTime();

    for (int field = ivx; field <= ietot; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid =i + (j + k*GridDimension[1])*GridDimension[0];
	    Prim[field][igrid] *= Prim[iden][igrid];
	    OldPrim[field][igrid] *= OldPrim[iden][igrid];
	  }
	}
      }
    }

    for (int field = 0; field < NEQ_MHD; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid =i + (j + k*GridDimension[1])*GridDimension[0];
	    Prim[field][igrid] = 0.5*(OldPrim[field][igrid] + Prim[field][igrid]);
	  }
	}
      }
    }

    for (int field = ivx; field <= ietot; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	    int igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	    Prim[field][igrid] /= Prim[iden][igrid];
	    OldPrim[field][igrid] /= OldPrim[iden][igrid];
	  }
	}
      }
    }

    return SUCCESS;

  } // if (UseCUDA)
#endif // ifdef ECUDA

  if (StellarWindFeedback)
    this->ReduceWindBoundary();

  float *dU[NEQ_MHD+NSpecies+NColor];
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
  }

  this->ReturnHydroRKPointers(Prim, true);  //##### added! because Hydro3D needs fractions for species

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

  return SUCCESS;

}
