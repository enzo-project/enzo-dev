/***********************************************************************
/
/  GRID CLASS (RUNGE-KUTTA FIRST STEP)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "EnzoTiming.h"
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


int grid::MHDRK2_1stStep(fluxes *SubgridFluxes[], 
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

  if (DualEnergyFormalism > 0) NEQ_MHD = 10;

  TIMER_START("MHDRK2");
#ifdef ECUDA
  if (UseCUDA) {
    this->CudaMHDRK2_1stStep(SubgridFluxes, NumberOfSubgrids, level, Exterior);
    TIMER_STOP("MHDRK2");
    return SUCCESS;
  }
#endif


  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
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

  float *Prim[NEQ_MHD+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false); 

  /* RK2 first step */


  float *dU[NEQ_MHD+NSpecies+NColor];

  int activesize = 1, i, dim;
  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++)
    dU[field] = new float[activesize];

  this->ReturnHydroRKPointers(Prim, true); //##### added! because Hydro3D needs fractions for species


  /* Compute dU */

  int fallback = 0;
  if (this->MHD3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 
		  0.5, fallback) == FAIL) {
    return FAIL;
  }

  /* Add source terms */

  this->MHDSourceTerms(dU);

  /* Update primitive variables */

  if (this->UpdateMHDPrim(dU, 1, 1) == FAIL) {
    fprintf(stderr, "Grid_MHDRK2_1stStep: Falling back to zero order at RK 1st step\n");
        // fall back to zero order scheme
    this->CopyOldBaryonFieldToBaryonField();
    for (int ns = NEQ_MHD; ns < NEQ_MHD+NSpecies+NColor; ns++) {
      // change species from density to mass fraction
      for (int n = 0; n < size; n++) {
	Prim[ns][n] /= Prim[iden][n];
      }
    }
    this->ZeroFluxes(SubgridFluxes, NumberOfSubgrids);
    fallback = 1;
    if (this->MHD3D(Prim, dU, dtFixed, SubgridFluxes, NumberOfSubgrids, 
                    0.5, fallback) == FAIL) {
      return FAIL;
    }
    this->MHDSourceTerms(dU);
    if (this->UpdateMHDPrim(dU, 1, 1) == FAIL) {
      fprintf(stderr, "Grid_MHDRK2_1stStep: Fallback failed, give up...\n");
      return FAIL;
    }
  }

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] dU[field];
  }

  TIMER_STOP("MHDRK2");

  return SUCCESS;

}
