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

int grid::MHDRK2_1stStep(fluxes *SubgridFluxes[], 
			 int NumberOfSubgrids, int level,
			 ExternalBoundary *Exterior)
  /*
    NumberOfSubgrids: the actual number of subgrids + 1
    SubgridFluxes[NumberOfSubgrids]
  */
{
  //  printf("NumberOfBaryonFields=%"ISYM"\n", NumberOfBaryonFields);
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

  if (DualEnergyFormalism > 0) NEQ_MHD = 10;

  float *Prim[NEQ_MHD+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false); 

  /* RK2 first step */


#ifdef ECUDA
  if (UseCUDA == 1) {
    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time2 = ReturnWallTime();
    if (MHDTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
			    dtdx, dtFixed, C_h, C_p, Theta_Limiter) == FAIL) {
      printf("RK1: MHDTimeUpdate_CUDA failed.\n");
      return FAIL;
    }
    return SUCCESS;
  }
#endif

  float *dU[NEQ_MHD+NSpecies+NColor];

  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);

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
    return FAIL;
  }

  for (int field = 0; field < NEQ_MHD+NSpecies+NColor; field++) {
    delete [] dU[field];
  }

  return SUCCESS;

}
