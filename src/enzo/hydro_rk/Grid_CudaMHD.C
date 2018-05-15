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
#include "CUDAUtil.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

extern dim3 CudaBlock, CudaGrid;

// Helper routine using pointer arithmetics for memory allocation
void grid::CudaMHDMalloc(void **p, size_t size)
{
  assert(size > 0);
  //CUDA_SAFE_CALL( cudaMalloc(p, size) );
  size_t offset = 128*(int)ceil((float)size/128);
  if (MHDData.GPUMemOffset + offset > MHDData.GPUMemSize) {
    printf("insufficient GPU memory: GPUMemOffset=%d, offset=%d, GPUMemSize=%d!\n",
           MHDData.GPUMemOffset, offset, MHDData.GPUMemSize);
    exit(1);
  }
  *p = (void*)((char*)MHDData.GPUMem + MHDData.GPUMemOffset);
  MHDData.GPUMemOffset += offset;
}

// allocate space for all GPU arrays
void grid::CudaMHDMallocGPUData()
{

  // Init parameters for GPU MHD solvers
  for (int i = 0; i < GridRank; i++) {
    MHDData.Dimension[i] = GridDimension[i];
    MHDData.StartIndex[i] = GridStartIndex[i];
    MHDData.EndIndex[i] = GridEndIndex[i];
  }
  size_t size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  CudaBlock.x = NBLOCK;
  CudaGrid.x = (size+CudaBlock.x-1)/CudaBlock.x;
  if (CudaGrid.x > 65535)
    CudaGrid.y = (CudaGrid.x+255)/256, CudaGrid.x = 256;

  const size_t sizebytes = size*sizeof(float);
  size_t sizebytes_align = 128*(int)ceil((float)sizebytes/128);

  // Compute the size of the total GPU memory
  MHDData.GPUMemOffset = 0;
  MHDData.GPUMemSize = 4*NEQ_MHD*sizebytes_align;
#ifdef DEDNER_SOURCE
  MHDData.GPUMemSize += 4*sizebytes_align;
#endif

  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    MHDData.GPUMemSize += 3*sizebytes_align;
  if (UseDrivingField) 
    MHDData.GPUMemSize += 3*sizebytes_align;
  if (MultiSpecies)
    MHDData.GPUMemSize += 4*NSpecies*sizebytes_align;

  size_t free, tot;
  CUDA_SAFE_CALL( cudaMemGetInfo(&free, &tot) );
  free -= 200*1024*1024;
  if (MHDData.GPUMemSize > free) {
    printf("Requested GPU memory size %f MB > available GPU memory %f MB\n",
           (float)MHDData.GPUMemSize/(1024.0*1024.0), 
           (float)free/(1024.0*1024.0));
    exit(1);
  }
  // allocate memory for all variables
  CUDA_SAFE_CALL( cudaMalloc(&(MHDData.GPUMem), MHDData.GPUMemSize) );

  // baryon 
  for (int i = 0; i < NEQ_MHD; i++) {
    CudaMHDMalloc((void**)&MHDData.Baryon[i], sizebytes);
    CudaMHDMalloc((void**)&MHDData.OldBaryon[i], sizebytes);
    CudaMHDMalloc((void**)&MHDData.Flux[i], sizebytes);
    CudaMHDMalloc((void**)&MHDData.dU[i], sizebytes);
  }

  // source terms
  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    for (int i = 0; i < GridRank; i++) {
      CudaMHDMalloc((void**)&MHDData.AccelerationField[i], sizebytes);
    }
  if (UseDrivingField) 
    for (int i = 0; i < GridRank; i++)
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
  cudaFree(MHDData.GPUMem);
  MHDData.GPUMemOffset = 0;
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
  // Set dU to zero if we are the first sweep
  if (dir == 1) {
    const int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
    const int sizebytes = size*sizeof(float); 
    for (int i = 0; i < NEQ_MHD; i++)
      cudaMemset(MHDData.dU[i], 0, sizebytes);
    for (int i = 0; i < NSpecies; i++)
      cudaMemset(MHDData.dUSpecies[i], 0, sizebytes);
  }
 
  // HLL-PLM solver to compute flux
  MHD_HLL_PLMGPU(MHDData, dir, CellWidth[0][0]);

  // Compute flux for species field
  if (NSpecies) 
    ComputeFluxSpeciesGPU(MHDData, dir);
 
  float a=1.0f, dadt=0.0f;
  if (ComovingCoordinates) 
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
  float dx = CellWidth[dir-1][0]*a;

  // Compute dU from flux difference
  ComputedUGPU(MHDData, dtFixed, dx, dir);
}

void grid::CudaMHDSourceTerm()
{
  float a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

  if (DualEnergyFormalism) 
    MHDDualEnergySourceGPU(MHDData, dtFixed, a,
                           CellWidth[0][0], CellWidth[1][0], CellWidth[2][0]);

  if (SelfGravity || ExternalGravity || UniformGravity || PointSourceGravity) 
    MHDGravitySourceGPU(MHDData, dtFixed);

  if (ComovingCoordinates) {
    float coef = -0.5f*dadt/a;
    MHDComovingSourceGPU(MHDData, dtFixed, coef);
  }
  
  if (UseDrivingField) 
    MHDDrivingSourceGPU(MHDData, dtFixed);
}

// RK2 update on GPU
void grid::CudaMHDUpdatePrim(int RKStep)
{  
  UpdateMHDPrimGPU(MHDData, RKStep, dtFixed);
}

void grid::CudaSolveMHDEquations(fluxes *SubgridFluxes[], int NumberOfSubgrids, int RKStep)
{
  const int size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  // Convert species from density to mass fraction
  if (NSpecies > 0) 
    Density2FractionGPU(MHDData);

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
  if (NSpecies > 0) 
    Fraction2DensityGPU(MHDData);
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
    int FluxId[10] = 
      {DensNum, Vel1Num, Vel2Num, Vel3Num, TENum,
       B1Num, B2Num, B3Num, PhiNum, GENum};
    for (int i = 0; i < NEQ_MHD; i++) {
      MHDSaveSubgridFluxGPU(SubgridFluxes[n]->LeftFluxes[FluxId[i]][dim],
                            SubgridFluxes[n]->RightFluxes[FluxId[i]][dim],
                            MHDData.Flux[i],
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
