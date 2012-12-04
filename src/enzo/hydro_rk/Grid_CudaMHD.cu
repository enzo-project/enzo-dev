/***********************************************************************
/
/  GRID CLASS (MHD SOLVER ON GPU)
/
/  written by: Peng Wang
/  date:       September, 2012
/  modified1:
/
/
************************************************************************/

 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "euler_sweep.h"
#include "fortran.def"

#include "CudaMHD.h"
#include "CudaMHD.cuh"
#include "../CUDAUtil.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

// Allocate a single memory space and use pointer arithmetics for all arrays
static void *GPUMem;
static size_t GPUMemSize;
static size_t GPUMemOffset;

static dim3 CudaBlock, CudaGrid;

// Helper routine using pointer arithmetics for memory allocation
void CudaMHDMalloc(void **p, size_t size)
{
  assert(size > 0);
  //CUDA_SAFE_CALL( cudaMalloc(p, size) );
  size_t offset = 128*(int)ceil((float)size/128);
  if (GPUMemOffset + offset > GPUMemSize) {
    printf("insufficient GPU memory!\n");
    exit(1);
  }
  *p = (void*)((char*)GPUMem + GPUMemOffset);
  GPUMemOffset += offset;
}

// allocate space for all GPU arrays
void grid::CudaMHDMallocGPUData()
{
  size_t size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  CudaBlock.x = NBLOCK;
  CudaGrid.x = (size+CudaBlock.x-1)/CudaBlock.x;
  if (CudaGrid.x > 65535)
    CudaGrid.y = (CudaGrid.x+255)/256, CudaGrid.x = 256;

  const size_t sizebytes = size*sizeof(float);
  size_t sizebytes_align = 128*(int)ceil((float)sizebytes/128);

  // Compute the size of the total GPU memory
  GPUMemSize = 36*sizebytes_align;
  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    GPUMemSize += 3*sizebytes_align;
  if (UseDrivingField) 
    GPUMemSize += 3*sizebytes_align;
  if (MultiSpecies)
    GPUMemSize += 4*NSpecies*sizebytes_align;

  size_t free, tot;
  CUDA_SAFE_CALL( cudaMemGetInfo(&free, &tot) );
  free -= 200*1024*1024;
  if (GPUMemSize > free) {
    printf("Requested GPU memory size %f MB > available GPU memory %f MB\n",
           (float)GPUMemSize/(1024.0*1024.0), (float)free/(1024.0*1024.0));
    exit(1);
  }
  CUDA_SAFE_CALL( cudaMalloc(&GPUMem, GPUMemSize) );
  
  // baryon 
  CudaMHDMalloc((void**)&MHDData.D  , sizebytes);
  CudaMHDMalloc((void**)&MHDData.V1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.V2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.V3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.TE , sizebytes);
  CudaMHDMalloc((void**)&MHDData.B1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.B2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.B3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.Phi, sizebytes);
  // old baryon
  CudaMHDMalloc((void**)&MHDData.OldD  , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldV1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldV2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldV3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldTE , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldB1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldB2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldB3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.OldPhi, sizebytes);
  // fluxes
  CudaMHDMalloc((void**)&MHDData.FluxD  , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxS1 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxS2 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxS3 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxTau, sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxB1 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxB2 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxB3 , sizebytes);  
  CudaMHDMalloc((void**)&MHDData.FluxPhi, sizebytes);  
  // dU
  CudaMHDMalloc((void**)&MHDData.dUD  , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUS1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUS2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUS3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUTau, sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUB1 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUB2 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUB3 , sizebytes);
  CudaMHDMalloc((void**)&MHDData.dUPhi, sizebytes);
  // source terms
  // CudaMHDMalloc((void**)&MHDData.divB, sizebytes);
  // CudaMHDMalloc((void**)&MHDData.gradPhi, 3*sizebytes);
  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    for (int i = 0; i < 3; i++)
      CudaMHDMalloc((void**)&MHDData.AccelerationField[i], sizebytes);
  if (UseDrivingField) 
    for (int i = 0; i < 3; i++)
      CudaMHDMalloc((void**)&MHDData.DrivingForce[i], sizebytes);
  if (MultiSpecies) {
    cudaMalloc(&MHDData.SpeciesArray, NSpecies*sizeof(float*));
    cudaMalloc(&MHDData.OldSpeciesArray, NSpecies*sizeof(float*));
    cudaMalloc(&MHDData.FluxSpeciesArray, NSpecies*sizeof(float*));
    cudaMalloc(&MHDData.dUSpeciesArray, NSpecies*sizeof(float*));
    for (int i = 0; i < NSpecies; i++) {
      CudaMHDMalloc((void**)&MHDData.Species[i], sizebytes);
      CudaMHDMalloc((void**)&MHDData.OldSpecies[i], sizebytes);
      CudaMHDMalloc((void**)&MHDData.FluxSpecies[i], sizebytes);
      CudaMHDMalloc((void**)&MHDData.dUSpecies[i], sizebytes);
    }
    cudaMemcpy(MHDData.SpeciesArray, MHDData.Species, NSpecies*sizeof(float*),
               cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.OldSpeciesArray, MHDData.OldSpecies, NSpecies*sizeof(float*),
               cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.FluxSpeciesArray, MHDData.FluxSpecies, NSpecies*sizeof(float*),
               cudaMemcpyHostToDevice);
    cudaMemcpy(MHDData.dUSpeciesArray, MHDData.dUSpecies, NSpecies*sizeof(float*),
               cudaMemcpyHostToDevice);
  }
}                       
 
// free space for all GPU arrays
void grid::CudaMHDFreeGPUData()
{
  cudaFree(GPUMem);
  GPUMemOffset = 0;
  if (MultiSpecies) {
    cudaFree(MHDData.SpeciesArray);
    cudaFree(MHDData.OldSpeciesArray);
    cudaFree(MHDData.FluxSpeciesArray);
    cudaFree(MHDData.dUSpeciesArray);
  }
}
 
// HLL-PLM solver on GPU
void grid::CudaMHDSweep(int dir)
{
  int idim = GridDimension[0];
  int jdim = GridDimension[1]; 
  int kdim = GridDimension[2];
  int i1 = GridStartIndex[0];
  int i2 = GridEndIndex[0];
  int j1 = GridStartIndex[1];
  int j2 = GridEndIndex[1];
  int k1 = GridStartIndex[2];
  int k2 = GridEndIndex[2];

  if (dir == 1) {
    const int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
    const int sizebytes = size*sizeof(float); 
    cudaMemset(MHDData.dUD, 0, sizebytes);
    cudaMemset(MHDData.dUS1, 0, sizebytes);
    cudaMemset(MHDData.dUS2, 0, sizebytes);
    cudaMemset(MHDData.dUS3, 0, sizebytes);
    cudaMemset(MHDData.dUTau, 0, sizebytes);
    cudaMemset(MHDData.dUB1, 0, sizebytes);
    cudaMemset(MHDData.dUB2, 0, sizebytes);
    cudaMemset(MHDData.dUB3, 0, sizebytes);
    cudaMemset(MHDData.dUPhi, 0, sizebytes);
    //cudaMemset(MHDData.divB, 0, sizebytes);
    for (int i = 0; i < NSpecies; i++)
      cudaMemset(MHDData.dUSpecies[i], 0, sizebytes);
  }
 
  // HLL-PLM solver to compute flux
  MHD_HLL_PLMKernel<<<CudaGrid, CudaBlock>>>(
    MHDData.FluxD, MHDData.FluxS1, MHDData.FluxS2, MHDData.FluxS3, MHDData.FluxTau, 
    MHDData.FluxB1, MHDData.FluxB2, MHDData.FluxB3, MHDData.FluxPhi,
    MHDData.D, MHDData.V1, MHDData.V2, MHDData.V3, MHDData.TE,
    MHDData.B1, MHDData.B2, MHDData.B3, MHDData.Phi,
    EOSType, C_h, Theta_Limiter, Gamma, NEQ_MHD,
    idim, jdim, kdim, i1, i2, j1, j2, k1, k2, dir);
  CUDA_SAFE_CALL( cudaGetLastError() );

  // Compute flux for species field
  if (NSpecies) {
    ComputeFluxSpeciesKernel<<<CudaGrid, CudaBlock>>>(
      MHDData.FluxSpeciesArray, MHDData.SpeciesArray, MHDData.FluxD,
      NSpecies, Theta_Limiter, 
      idim, jdim, kdim, i1, i2, j1, j2, k1, k2, dir);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }
 
  float a=1.0f, dadt=0.0f;
  if (ComovingCoordinates) 
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
  float dx = CellWidth[dir-1][0]*a;

  // Compute dU from flux difference
  ComputedUKernel<<<CudaGrid, CudaBlock>>>(
    MHDData.dUD, MHDData.dUS1, MHDData.dUS2, MHDData.dUS3, MHDData.dUTau,
    MHDData.dUB1, MHDData.dUB2, MHDData.dUB3, MHDData.dUPhi, MHDData.divB, 
    MHDData.gradPhi, MHDData.dUSpeciesArray,
    MHDData.FluxD, MHDData.FluxS1, MHDData.FluxS2, MHDData.FluxS3, 
    MHDData.FluxTau, MHDData.FluxB1, MHDData.FluxB2, MHDData.FluxB3, 
    MHDData.FluxPhi, MHDData.FluxSpeciesArray,
    NSpecies, dtFixed, dx, C_h, idim, jdim, kdim, 
    i1, i2, j1, j2, k1, k2, dir);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

void grid::CudaMHDSourceTerm()
{
  float a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) {  
    MHDGravitySourceKernel<<<CudaGrid, CudaBlock>>>
      (MHDData.dUS1, MHDData.dUS2, MHDData.dUS3, MHDData.dUTau,
       MHDData.D, MHDData.V1, MHDData.V2, MHDData.V3, dtFixed,
       MHDData.AccelerationField[0], 
       MHDData.AccelerationField[1],
       MHDData.AccelerationField[2],
       GridDimension[0], GridDimension[1], GridDimension[2],
       GridStartIndex[0], GridEndIndex[0],
       GridStartIndex[1], GridEndIndex[1],
       GridStartIndex[2], GridEndIndex[2]);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }

  if (ComovingCoordinates) {
    float coef = -0.5f*dadt/a;
    MHDComovingSourceKernel<<<CudaGrid, CudaBlock>>>
      (MHDData.dUTau, MHDData.dUB1, MHDData.dUB2, MHDData.dUB3, MHDData.dUPhi,
       MHDData.B1, MHDData.B2, MHDData.B3, dtFixed, coef,
       GridDimension[0], GridDimension[1], GridDimension[2],
       GridStartIndex[0], GridEndIndex[0],
       GridStartIndex[1], GridEndIndex[1],
       GridStartIndex[2], GridEndIndex[2]);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }
  
  if (UseDrivingField) {
    MHDDrivingSourceKernel<<<CudaGrid, CudaBlock>>>
      (MHDData.dUS1, MHDData.dUS2, MHDData.dUS3, MHDData.dUTau,
       MHDData.D, MHDData.V1, MHDData.V2, MHDData.V3, MHDData.TE,
       MHDData.B1, MHDData.B2, MHDData.B3,
       MHDData.DrivingForce[0], MHDData.DrivingForce[1], MHDData.DrivingForce[2],
       dtFixed, DrivingEfficiency,
       GridDimension[0], GridDimension[1], GridDimension[2],
       GridStartIndex[0], GridEndIndex[0],
       GridStartIndex[1], GridEndIndex[1],
       GridStartIndex[2], GridEndIndex[2]);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }       
}

// RK2 update on GPU
void grid::CudaMHDUpdatePrim(int RKStep)
{  
  UpdateMHDPrimKernel<<<CudaGrid, CudaBlock>>>(
    MHDData.D, MHDData.V1, MHDData.V2, MHDData.V3, MHDData.TE, 
    MHDData.B1, MHDData.B2, MHDData.B3, MHDData.Phi, MHDData.SpeciesArray,
    MHDData.OldD, MHDData.OldV1, MHDData.OldV2, MHDData.OldV3, MHDData.OldTE,
    MHDData.OldB1, MHDData.OldB2, MHDData.OldB3, MHDData.OldPhi, MHDData.OldSpeciesArray,
    MHDData.dUD, MHDData.dUS1, MHDData.dUS2, MHDData.dUS3, MHDData.dUTau,
    MHDData.dUB1, MHDData.dUB2, MHDData.dUB3, MHDData.dUPhi, MHDData.dUSpeciesArray,
    NSpecies, RKStep, dtFixed, C_h, C_p, 
    GridDimension[0], GridDimension[1], GridDimension[2],
    GridStartIndex[0], GridEndIndex[0],
    GridStartIndex[1], GridEndIndex[1],
    GridStartIndex[2], GridEndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

void grid::CudaSolveMHDEquations(fluxes *SubgridFluxes[], int NumberOfSubgrids, int RKStep)
{
  const int size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  // Convert species from density to mass fraction
  if (NSpecies > 0) {
    Density2Fraction<<<CudaGrid, CudaBlock>>>(MHDData.SpeciesArray, MHDData.D, NSpecies, size);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }

  this->CudaMHDSweep(1);

  if (FluxCorrection) 
    this->CudaMHDSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 1);
  
  this->CudaMHDSweep(2);

  if (FluxCorrection) 
    this->CudaMHDSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 2); 

  this->CudaMHDSweep(3);

  if (FluxCorrection) 
    this->CudaMHDSaveSubgridFluxes(SubgridFluxes, NumberOfSubgrids, 3); 

  this->CudaMHDSourceTerm();

  this->CudaMHDUpdatePrim(RKStep);

  // Convert species from mass fraction back to density
  if (NSpecies > 0) {
    Fraction2Density<<<CudaGrid, CudaBlock>>>(MHDData.SpeciesArray, MHDData.D, NSpecies, size);
    CUDA_SAFE_CALL( cudaGetLastError() );
  }
}

int MHDSaveSubgridFluxCUDA(float *LeftFlux, float *RightFlux, float *Flux,
                           float *Tmp1, float *Tmp2, float *Tmp3, float *Tmp4, 
                           float dtdx,
                           const int dimx, const int dimy, const int dimz,
                           const int fistart, const int fiend,
                           const int fjstart, const int fjend,
                           const int lface, const int rface, int dir)
{
  int nf = (fiend - fistart + 1)*(fjend - fjstart + 1);
  
  MHDSaveSubgridFluxKernel<<<(nf+127)/128, 128>>>
    (Tmp1, Tmp2, Flux, dtdx, dimx, dimy, dimz,
     fistart, fiend, fjstart, fjend, lface, rface, dir);
  cudaMemcpy(Tmp3, Tmp1, nf*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(Tmp4, Tmp2, nf*sizeof(float), cudaMemcpyDeviceToHost);  
  for (int i = 0; i < nf; i++) {
    LeftFlux[i] += Tmp3[i];
    RightFlux[i] += Tmp4[i];
  }
  return SUCCESS;
}                     

void grid::CudaMHDSaveSubgridFluxes(fluxes *SubgridFluxes[],
                                    int NumberOfSubgrids,
                                    int dir)
{
  float a=1.0f, dadt=0.0f;
  if (ComovingCoordinates) 
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
  float dtdx = dtFixed/(CellWidth[dir-1][0]*a);

  int dim, dimi, dimj;
  if (dir == 1)
    dim = 0, dimi = 1, dimj = 2;
  else if (dir == 2)
    dim = 1, dimi = 0, dimj = 2;
  else if (dir == 3)
    dim = 2, dimi = 0, dimj = 1;

  Elong_int GridGlobalStart[MAX_DIMENSION];
  for (int dim = 0; dim < GridRank; dim++)
    GridGlobalStart[dim] = 
      nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
	GridStartIndex[dim];

  for (int n = 0; n < NumberOfSubgrids; n++) {
    const int fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dimi] -
      GridGlobalStart[dimi];
    const int fiend   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][dimi] -
      GridGlobalStart[dimi];
    const int fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dimj] - 
      GridGlobalStart[dimj];
    const int fjend   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][dimj] -
      GridGlobalStart[dimj];
    const int lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] - 
      GridGlobalStart[dim];
    const int rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
      GridGlobalStart[dim] + 1;
    int nfi = fiend - fistart + 1;
    int nfj = fjend - fjstart + 1;
    int nf = nfi*nfj;
    float *LeftFlux, *RightFlux;
    cudaMalloc(&LeftFlux, nf*sizeof(float));
    cudaMalloc(&RightFlux, nf*sizeof(float));
    float *LeftTmp = (float*)malloc(nf*sizeof(float));
    float *RightTmp = (float*)malloc(nf*sizeof(float));
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                     Vel3Num, TENum, B1Num, B2Num, B3Num, 
                                     PhiNum);
    int FluxId[9] = 
      {DensNum, Vel1Num, Vel2Num, Vel3Num, TENum,
       B1Num, B2Num, B3Num, PhiNum};
    float *Flux3D[9] = 
      {MHDData.FluxD, MHDData.FluxS1, MHDData.FluxS2, MHDData.FluxS3, 
       MHDData.FluxTau, MHDData.B1, MHDData.B2, MHDData.B3, MHDData.Phi};
    for (int i = 0; i < 9; i++) {
      MHDSaveSubgridFluxCUDA(SubgridFluxes[n]->LeftFluxes[FluxId[i]][dim],
                             SubgridFluxes[n]->RightFluxes[FluxId[i]][dim],
                             Flux3D[i],
                             LeftFlux, RightFlux, LeftTmp, RightTmp, dtdx,
                             GridDimension[0],
                             GridDimension[1],
                             GridDimension[2],
                             fistart, fiend, fjstart, fjend, lface, rface, dir);
    }
    cudaFree(LeftFlux);
    cudaFree(RightFlux);
    free(LeftTmp);
    free(RightTmp);
  }
}
