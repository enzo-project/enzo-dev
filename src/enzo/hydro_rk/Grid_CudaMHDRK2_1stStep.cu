/***********************************************************************
/
/  GRID CLASS (RUNGE-KUTTA FIRST STEP ON GPU)
/
/  written by: Peng Wang
/  date:       September, 2012
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "CUDAUtil.h"


int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::CudaMHDRK2_1stStep(fluxes *SubgridFluxes[], 
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

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

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

  /* RK2 first step */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                   Vel3Num, TENum, B1Num, B2Num, B3Num, 
                                   PhiNum);
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies)
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  const size_t sizebytes = size*sizeof(float);
  //
  // Allocate memory on GPU
  //
  this->CudaMHDMallocGPUData();

  //
  // Copy data from CPU to GPU
  //
  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    for (int i = 0; i < GridRank; i++)
      cudaMemcpy(MHDData.AccelerationField[i], AccelerationField[i], sizebytes, cudaMemcpyHostToDevice);
  if (UseDrivingField) {
    int Drive1Num, Drive2Num, Drive3Num;
    if (IdentifyDrivingFields(Drive1Num, Drive2Num, Drive3Num) == FAIL) {
      printf("grid::CudaMHDRK2_1stStep: canot identify driving fields.\n");
      return FAIL;
    }
    cudaMemcpy(MHDData.DrivingForce[0], BaryonField[Drive1Num], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.DrivingForce[1], BaryonField[Drive2Num], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.DrivingForce[2], BaryonField[Drive3Num], sizebytes, cudaMemcpyHostToDevice);
  }   
  cudaMemcpy(MHDData.D  , BaryonField[DensNum], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.V1 , BaryonField[Vel1Num], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.V2 , BaryonField[Vel2Num], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.V3 , BaryonField[Vel3Num], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.TE , BaryonField[TENum  ], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.B1 , BaryonField[B1Num  ], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.B2 , BaryonField[B2Num  ], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.B3 , BaryonField[B3Num  ], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(MHDData.Phi, BaryonField[PhiNum ], sizebytes, cudaMemcpyHostToDevice);

  // copy to old baryon
  cudaMemcpy(MHDData.OldD  , MHDData.D  , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldV1 , MHDData.V1 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldV2 , MHDData.V2 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldV3 , MHDData.V3 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldTE , MHDData.TE , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldB1 , MHDData.B1 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldB2 , MHDData.B2 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldB3 , MHDData.B3 , sizebytes, cudaMemcpyDeviceToDevice);
  cudaMemcpy(MHDData.OldPhi, MHDData.Phi, sizebytes, cudaMemcpyDeviceToDevice);

  if (MultiSpecies) {
    cudaMemcpy(MHDData.Species[0], BaryonField[HINum], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.Species[1], BaryonField[HIINum], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.Species[2], BaryonField[HeINum], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.Species[3], BaryonField[HeIINum], sizebytes, cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.Species[4], BaryonField[HeIIINum], sizebytes, cudaMemcpyHostToDevice);
    if (MultiSpecies > 1) {
      cudaMemcpy(MHDData.Species[5], BaryonField[HMNum], sizebytes, cudaMemcpyHostToDevice);
      cudaMemcpy(MHDData.Species[6], BaryonField[H2INum], sizebytes, cudaMemcpyHostToDevice);
      cudaMemcpy(MHDData.Species[7], BaryonField[H2IINum], sizebytes, cudaMemcpyHostToDevice);
    }
    if (MultiSpecies > 2) {
      cudaMemcpy(MHDData.Species[8], BaryonField[DINum], sizebytes, cudaMemcpyHostToDevice);
      cudaMemcpy(MHDData.Species[9], BaryonField[DIINum], sizebytes, cudaMemcpyHostToDevice);
      cudaMemcpy(MHDData.Species[10], BaryonField[HDINum], sizebytes, cudaMemcpyHostToDevice);
    }
    // copy to old species
    for (int i = 0; i < NSpecies; i++)
      cudaMemcpy(MHDData.OldSpecies[i], MHDData.Species[i], sizebytes,
                 cudaMemcpyDeviceToDevice);
  }
  CUDA_SAFE_CALL( cudaGetLastError() );
  //
  // Solve MHD equations on GPU
  //
  this->CudaSolveMHDEquations(SubgridFluxes, NumberOfSubgrids, 1);

  //
  // Copy results from CPU to GPU                                                    
  //
  cudaMemcpy(BaryonField[DensNum], MHDData.D  , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Vel1Num], MHDData.V1 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Vel2Num], MHDData.V2 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Vel3Num], MHDData.V3 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[TENum  ], MHDData.TE , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[B1Num  ], MHDData.B1 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[B2Num  ], MHDData.B2 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[B3Num  ], MHDData.B3 , sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[PhiNum ], MHDData.Phi, sizebytes, cudaMemcpyDeviceToHost);
  if (MultiSpecies) {
    cudaMemcpy(BaryonField[HINum   ], MHDData.Species[0], sizebytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(BaryonField[HIINum  ], MHDData.Species[1], sizebytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(BaryonField[HeINum  ], MHDData.Species[2], sizebytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(BaryonField[HeIINum ], MHDData.Species[3], sizebytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(BaryonField[HeIIINum], MHDData.Species[4], sizebytes, cudaMemcpyDeviceToHost);
    if (MultiSpecies > 1) {
      cudaMemcpy(BaryonField[HMNum  ], MHDData.Species[5], sizebytes, cudaMemcpyDeviceToHost);
      cudaMemcpy(BaryonField[H2INum ], MHDData.Species[6], sizebytes, cudaMemcpyDeviceToHost);
      cudaMemcpy(BaryonField[H2IINum], MHDData.Species[7], sizebytes, cudaMemcpyDeviceToHost);
    }
    if (MultiSpecies > 2) {
      cudaMemcpy(BaryonField[DINum ], MHDData.Species[8], sizebytes, cudaMemcpyDeviceToHost);
      cudaMemcpy(BaryonField[DIINum], MHDData.Species[9], sizebytes, cudaMemcpyDeviceToHost);
      cudaMemcpy(BaryonField[HDINum], MHDData.Species[10],sizebytes, cudaMemcpyDeviceToHost);
    }
  }
  if (NSpecies > 0)
    this->UpdateElectronDensity();

  return SUCCESS;
}

