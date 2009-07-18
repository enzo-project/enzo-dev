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

int grid::MHDRK2_1stStep(int CycleNumber, fluxes *SubgridFluxes[], 
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

  // Dual Energy Formalism not fully implemented with MHD yet ... 
  
  if (DualEnergyFormalism > 0) NEQ_MHD = 10.;

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
  
  /* RK2 first step */

  /* Compute dU */

#ifdef ECUDA
  if (UseCUDA == 1) {
    FLOAT dtdx = dtFixed/CellWidth[0][0];
    double time2 = ReturnWallTime();
    if (MHDTimeUpdate_CUDA(Prim, GridDimension, GridStartIndex, GridEndIndex, GridRank,
			    dtdx, dtFixed, C_h, C_p) == FAIL) {
      printf("RK1: MHDTimeUpdate_CUDA failed.\n");
      return FAIL;
    }
    //return FAIL;
    //    PerformanceTimers[1] += ReturnWallTime() - time1;
    return SUCCESS;
  }
#endif


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

  //  PerformanceTimers[1] += ReturnWallTime() - time1;

  return SUCCESS;

}
