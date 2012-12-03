/***********************************************************************
/
/  GRID CLASS (FIRST STEP OF RUNGE-KUTTA INTEGRATION)
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
int HydroTimeUpdate_CUDA(float **Prim, int GridDimension[], 
			 int GridStartIndex[], int GridEndIndex[], int GridRank,
			 float dtdx, float dt);

int grid::RungeKutta2_1stStep(fluxes *SubgridFluxes[], 
			      int NumberOfSubgrids, int level,
			      ExternalBoundary *Exterior)  {
  /*
    NumberOfSubgrids: the actual number of subgrids + 1
    SubgridFluxes[NumberOfSubgrids]
  */
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (NumberOfBaryonFields == 0) {
    return SUCCESS;
  }

  double time1 = ReturnWallTime();
  int igrid;
  /* allocate space for fluxes */
  int fluxsize;
  for (int subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    for (int flux = 0; flux < GridRank; flux++)  {
      
      fluxsize = 1;
      for (int j = 0; j < GridRank; j++) {
	fluxsize *= SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[flux][j] -
	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][j] + 1;
      }
      
      for (int j = GridRank; j < 3; j++) {
	SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][j] = 0;
	SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[flux][j] = 0;
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[flux][j] = 0;
	SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[flux][j] = 0;
      }
       
      for (int field = 0; field < NumberOfBaryonFields; field++) {
	if (SubgridFluxes[subgrid]->LeftFluxes[field][flux] == NULL) {
	  SubgridFluxes[subgrid]->LeftFluxes[field][flux]  = new float[fluxsize];
	}
	if (SubgridFluxes[subgrid]->RightFluxes[field][flux] == NULL)
	  SubgridFluxes[subgrid]->RightFluxes[field][flux] = new float[fluxsize];
	for (int n = 0; n < fluxsize; n++) {
	  SubgridFluxes[subgrid]->LeftFluxes[field][flux][n] = 0.0;
	  SubgridFluxes[subgrid]->RightFluxes[field][flux][n] = 0.0;
	}
      }
      
      for (int field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][flux] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][flux] = NULL;
      }
      
    }  // next flux
    
    for (int flux = GridRank; flux < 3; flux++) {
      for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][flux] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][flux] = NULL;
      }
    }
    
  } // end of loop over subgrids


  float *Prim[NEQ_HYDRO+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false);  

  // RK2 first step
#ifdef ECUDA 
  if (UseCUDA == 1) {
    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time2 = ReturnWallTime();
    if (HydroTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
		  	      dtdx, dtFixed) == FAIL) {
      printf("RK1: HydroTimeUpdate_CUDA failed.\n");
      return FAIL;
    }
    return SUCCESS;
  }
#endif

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

  if (this->UpdatePrim(dU, 1.0, 1.0) == FAIL) {
    printf("Falling back to zero order at RK 1st step\n");
    // fall back to zero order scheme
    this->CopyOldBaryonFieldToBaryonField();
    for (int ns = NEQ_HYDRO; ns < NEQ_HYDRO+NSpecies+NColor; ns++) {
      // change species from density to mass fraction
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
    if (this->UpdatePrim(dU, 1.0, 1.0) == FAIL) {
      printf("Fallback failed, give up...\n");
      return FAIL;
    }
  }

  for (int field = 0; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    delete [] dU[field];
  }
  //  PerformanceTimers[1] += ReturnWallTime() - time1;
  return SUCCESS;

}
