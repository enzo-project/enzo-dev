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
int HydroTimeUpdate_CUDA(float **Prim, int GridDimension[], int GridStartIndex[], int GridEndIndex[], int GridRank,
		          float dtdx, float dt);


int grid::RungeKutta2_2ndStep(fluxes *SubgridFluxes[], 
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

  float *Prim[NEQ_HYDRO+NSpecies+NColor];
  float *OldPrim[NEQ_HYDRO+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false); 
  this->ReturnOldHydroRKPointers(OldPrim, false); 

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
	    Prim[field][igrid] *= Prim[iden][igrid];
	    OldPrim[field][igrid] *= OldPrim[iden][igrid];
	  }
	}
      }
    }

    for (int field = 0; field < NEQ_HYDRO; field++) {
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

  } // if UseCUDA == 1
#endif /* ECUDA */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);

  float *dU[NEQ_HYDRO+NSpecies+NColor];
  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    dU[field] = new float[activesize];
    for (int i = 0; i < activesize; i++) {
      dU[field][i] = 0.0;
    }
  }

  this->ReturnHydroRKPointers(Prim, true);  //##### added! because Hydro3D needs fractions for species

  // compute dU
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
    for (int ns = NEQ_HYDRO; ns < NEQ_HYDRO+NSpecies+NColor; ns++) {
      for (int n = 0; n < size; n++) {
	Prim[ns][n] /= Prim[iden][n];
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
  
  /* If we're supposed to be outputting on Density, we need to update
  the current maximum value of that Density. */
  
  if(OutputOnDensity == 1){
    int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
    for(int i = 0; i < size; i++)
      CurrentMaximumDensity = max(BaryonField[DensNum][i], CurrentMaximumDensity);
  }

  return SUCCESS;

}
