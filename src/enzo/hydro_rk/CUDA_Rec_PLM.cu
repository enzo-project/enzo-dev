/***********************************************************************
/
/  PLM RECONSTRUCTION
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/
************************************************************************/

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "EOS.h"
//#include <math.h>
//#include <cuda.h>
//#include "ReconstructionRoutines.h"
#define GHOST_ZONES 3

#define FABS(A) (A*sign(A))
__device__  float minmod(float a, float b, float c)
{
  return 0.25*(sign(a)+sign(b))*FABS((sign(a)+sign(c)))*min(min(FABS(a), FABS(b)), FABS(c));
}


__device__  void cuda_plm_l(float &vm1, float &v, float &vp1, float &vl_plm, float Theta_Limiter)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

__device__  void cuda_plm_r(float &vm1, float &v, float &vp1, float &vr_plm, float Theta_Limiter)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vr_plm = v - 0.5*dv;
}

__device__  void cuda_plm_point(float &vm1, float &v, float &vp1, float &vl_plm, float Theta_Limiter)
{

  float dv_l, dv_r, dv_m, dv;
  
  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);
  
  dv = minmod(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

__global__ void cuda_plm_kernel(float **prim, float **priml, float **primr, int ActiveSize, int field, float Theta_Limiter)
{

  int iprim;

  int i = blockIdx.x * blockDim.x + threadIdx.x; 

   if ( i < ActiveSize+1) {
     iprim = i + GHOST_ZONES - 1;
     cuda_plm_point(prim[field][iprim-1], prim[field][iprim  ], prim[field][iprim+1],priml[field][i], Theta_Limiter);
     cuda_plm_point(prim[field][iprim+2], prim[field][iprim+1], prim[field][iprim], primr[field][i], Theta_Limiter);
   }
}


__host__ void cuda_plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq, float Theta_Limiter)
{
  int blocksize = 16;
  dim3 dimBlock(blocksize);
  dim3 dimGrid( ceil( ActiveSize / (float)blocksize)  );
  int Nzones = ActiveSize+2*GHOST_ZONES;
  float *prim_d[Nzones],*priml_d[Nzones],*primr_d[Nzones]; 
  cudaMalloc( (void**)&prim_d, Nzones*sizeof(float));
  cudaMalloc( (void**)&priml_d, Nzones*sizeof(float));
  cudaMalloc( (void**)&priml_d, Nzones*sizeof(float));
  for (int field = 0; field < Neq; field++) {
    (cudaMemcpy( prim_d,   prim[field], Nzones*sizeof(float),cudaMemcpyHostToDevice));
    (cudaMemcpy( priml_d, priml[field], Nzones*sizeof(float),cudaMemcpyHostToDevice ));
    (cudaMemcpy( primr_d, primr[field], Nzones*sizeof(float),cudaMemcpyHostToDevice ));
    cuda_plm_kernel<<< dimGrid,dimBlock >>>(prim_d, priml_d, primr_d, ActiveSize, field, Theta_Limiter);
    // only copy left and right values back 
    cudaMemcpy( priml_d, priml[field], Nzones*sizeof(float),cudaMemcpyDeviceToHost );
    cudaMemcpy( primr_d, primr[field], Nzones*sizeof(float),cudaMemcpyDeviceToHost );
  }
  cudaFree(  prim_d );
  cudaFree( priml_d );
  cudaFree( primr_d );

}
__device__ int cuda_plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize, float Theta_Limiter, int NSpecies)
  /* is : starting index of species field in prim */
{

  int iprim;
  for (int n = 0; n < ActiveSize+1; n++) {
    iprim = n + GHOST_ZONES - 1;    
    if (flux0[n] >= 0) {
      for (int field = 0; field < NSpecies; field++) {
	cuda_plm_l(prim[field+is][iprim-1], prim[field+is][iprim], prim[field+is][iprim+1],
	      species[field][n], Theta_Limiter);
      }
    } else {
      for (int field = 0; field < NSpecies; field++) {
	cuda_plm_r(prim[field+is][iprim], prim[field+is][iprim+1], prim[field+is][iprim+2],
	      species[field][n], Theta_Limiter);
      }
    }

    /* renormalize species field */
    float sum = 0;
    for (int field = 0; field < NSpecies; field++) {
      sum += species[field][n];
    }
    for (int field = 0; field < NSpecies; field++) {
      species[field][n] /= sum;
    }
  }

  return SUCCESS;
  
}

__device__ int cuda_plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize, float Theta_Limiter, int NSpecies, int NColor)
  /* is : starting index of *species* field in prim */
{

  int iprim;
  for (int n = 0; n < ActiveSize+1; n++) {
    iprim = n + GHOST_ZONES - 1;        
    if (flux0[n] >= 0) {
      for (int field = is+NSpecies; field < is+NSpecies+NColor; field++) {
	cuda_plm_l(prim[field][iprim-1], prim[field][iprim], prim[field][iprim+1],
		   color[field-is-NSpecies][n], Theta_Limiter);
      }
    } else {
      for (int field = is+NSpecies; field < is+NSpecies+NColor; field++) {
	cuda_plm_r(prim[field][iprim], prim[field][iprim+1], prim[field][iprim+2],
	      color[field-is-NSpecies][n], Theta_Limiter);
      }
    }

  }

  return SUCCESS;
  
}
