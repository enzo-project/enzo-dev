/***********************************************************************
/
/  SAVE SUBGRID FLUX ON GPU
/
/  written by: Peng Wang 
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/         Routines for saving subgrid flux on GPU
/
************************************************************************/

//#include "macros_and_parameters.h"
#include "CUDAUtil.h"
 
__global__ void SaveSubgridFluxKernel(float *LeftFluxes,
                                      float *RightFluxes,
                                      const float* __restrict Flux,
                                      int dimx, int dimy, int dimz,
                                      int fistart, int fiend,
                                      int fjstart, int fjend,
                                      int lface, int rface, int dir)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  const int nfi = fiend - fistart + 1;
  const int nfj = fjend - fjstart + 1;
  int i, j;
  if (tid < nfi*nfj) {
    i = (tid % nfi) + fistart;
    j = tid / nfi + fjstart;
    int offset = (i-fistart) + (j-fjstart)*nfi;
    int lindex, rindex;
    if (dir == 0) {
      lindex = i*dimx + j*dimx*dimy + lface;
      rindex = i*dimx + j*dimx*dimy + rface;
    } else if (dir == 1) {
      lindex = j*dimy + i*dimy*dimz + lface;
      rindex = j*dimy + i*dimy*dimz + rface;
    } else if (dir == 2) {
      lindex = i*dimz + j*dimz*dimx + lface;
      rindex = i*dimz + j*dimz*dimx + rface;
    }
    LeftFluxes[offset] = Flux[lindex];
    RightFluxes[offset] = Flux[rindex];
  }
}

void SaveSubgridFluxCUDA(float *LeftFlux, float *RightFlux, float *Flux,
                         float *Tmp1, float *Tmp2,
                         const int dimx, const int dimy, const int dimz,
                         const int fistart, const int fiend,
                         const int fjstart, const int fjend,
                         const int lface, const int rface, const int dir)
{
  int nf = (fiend - fistart + 1)*(fjend - fjstart + 1);

  SaveSubgridFluxKernel<<<(nf+127)/128, 128>>>
    (Tmp1, Tmp2, Flux, dimx, dimy, dimz, 
     fistart, fiend, fjstart, fjend, lface, rface, dir);
  cudaMemcpy(LeftFlux, Tmp1, nf*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(RightFlux, Tmp2, nf*sizeof(float), cudaMemcpyDeviceToHost);  
  //  CUDA_SAFE_CALL( cudaGetLastError() );
}                     

