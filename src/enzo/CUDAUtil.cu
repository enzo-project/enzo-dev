/***********************************************************************
/
/  GPU ROUTINES
/
/  written by: Peng Wang 
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/         Init GPU.
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "CUDAUtil.h"

extern "C"
int InitGPU(int ProcessRank)
{
  int NumGPU;
  CUDA_SAFE_CALL( cudaGetDeviceCount(&NumGPU) );
  if (NumGPU < 1) {
    printf("No NVIDIA GPU found. Cannot use CUDA version.\n");
    exit(1);
  }
  CUDA_SAFE_CALL( cudaSetDevice( ProcessRank % NumGPU ) );
  int *a;
  cudaMalloc(&a, sizeof(int));
  return SUCCESS;
}
