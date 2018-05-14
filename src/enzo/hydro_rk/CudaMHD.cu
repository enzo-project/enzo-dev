/***********************************************************************
/
/  MHD SOLVER ON GPU
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
#include "CUDAUtil.h"
#include "CudaMHD.h"

#define sign(A)  ((A) >  0  ?  1  : -1 )

dim3 CudaBlock, CudaGrid;

/***********************************************************
/
/   Function EOS: compute equation of state
/
/   Input: 
/     eostype: 
/       0: ideal gas
/     mode:  
/       1: given p and rho, calculate others.
/       2: given rho and e, calculate others.
/
************************************************************/

__forceinline__ __device__ 
void EOSDevice(float &p, float &rho, float &e, float &h, float &cs, 
               const float &Gamma, const int &eostype, const int &mode)
{
  float poverrho;
  
  if (eostype == 0) {
    
    if (mode == 1) {
      poverrho = p / rho;
      e = p / rho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e;
      poverrho = p / rho;
    }

    h = e + poverrho;
    cs = sqrt(Gamma*poverrho);
  }

  // if (eostype == 3) { // straight isothermal
  //   cs = EOSSoundSpeed;
  //   p = rho*cs*cs;
  //   e = p / ((cGamma-1.0)*rho);
  // }

}

// return the maximum of a, b, c
__forceinline__ __device__ 
float Max3Device(const float &a, const float &b, const float &c)  
{
  if (a > b) {
    if (a > c)
      return a;
    else 
      return c;
  } else {
    if (b > c)
      return b;
    else
      return c;
  }
}


// return the minimum of a, b, c
__forceinline__ __device__ 
float Min3Device(const float &a, const float &b, const float &c)
{
  if (a<b) {
    if (c<a)
      return c;
    else 
      return a;
  } else {
    if (c<b)
      return c;
    else 
      return b;
  }
}

#define FABS(A) ((A)*sign(A))

// minmod limiter function
__forceinline__ __device__ 
float MinmodDevice(const float &a, const float &b, const float &c)
{
  return 0.25f*(sign(a)+sign(b))*
    FABS((sign(a)+sign(c)))*Min3Device(FABS(a), FABS(b), FABS(c));
}


// do PLM reconstruction for a point
__forceinline__ __device__ 
void PLMPointDevice(float &vl_plm, 
                    const float &vm1, const float &v, const float &vp1, 
                    const float &ThetaLimiter)
{
  float dv_l = (v-vm1) * ThetaLimiter;
  float dv_r = (vp1-v) * ThetaLimiter;
  float dv_m = 0.5f*(vp1-vm1);
  
  float dv = MinmodDevice(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5f*dv;
  //vl_plm = ThetaLimiter;
}

__forceinline__ __device__ 
void PLMLeftDevice(float &vl_plm, 
                   const float &vm1, const float &v, const float &vp1,
                   const float &ThetaLimiter)
{
  float dv_l = (v-vm1) * ThetaLimiter;
  float dv_r = (vp1-v) * ThetaLimiter;
  float dv_m = 0.5f*(vp1-vm1);
  
  float dv = MinmodDevice(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5f*dv;
}

__forceinline__ __device__ 
void PLMRightDevice(float &vr_plm, 
                    const float &vm1, const float &v, const float &vp1,
                    const float &ThetaLimiter)
{
  float dv_l = (v-vm1) * ThetaLimiter;
  float dv_r = (vp1-v) * ThetaLimiter;
  float dv_m = 0.5f*(vp1-vm1);
  
  float dv = MinmodDevice(dv_l, dv_r, dv_m);

  vr_plm = v - 0.5f*dv;
}

// HLL Riemann solver
__forceinline__ __device__ 
void HLLDevice(float &UD, float &US1, float &US2, float &US3, 
               float &UTau, float &UB1, float &UB2, float &UB3, 
               float &UPhi, float &UGE,
               float &FD, float &FS1, float &FS2, float &FS3, 
               float &FTau, float &FB1, float &FB2, float &FB3, 
               float &FPhi, float &FGE,
               float &lp, float &lm,
               float D, float Ein, float V1, float V2, float V3, 
               float B1, float B2, float B3, float Phi,
               int EOSType, float Gamma, int DualEnergyFormalism, int idx3d)
{
  float BSq = B1*B1 + B2*B2 + B3*B3;
  float BV = B1*V1 + B2*V2 + B3*V3;
  float VSq = V1*V1 + V2*V2 + V3*V3;
  float TE = Ein + 0.5f*VSq + 0.5f*BSq/D;
  float p, h, cs;
  EOSDevice(p, D, Ein, h, cs, Gamma, EOSType, 2);
  if (EOSType > 0) {
    p = Ein;
    cs = sqrt(p/D);
  }
  float cs2 = cs*cs;

  UD  = D;
  US1 = D * V1;
  US2 = D * V2;
  US3 = D * V3;
  UTau = D * TE;
  UB1 = B1;
  UB2 = B2;
  UB3 = B3;
  UPhi = Phi;
  if (DualEnergyFormalism)
    UGE = D * Ein;

  FD = D * V1;
  FS1 = US1 * V1 + p + 0.5f*BSq - B1*B1;
  FS2 = US2 * V1 - B1*B2;
  FS3 = US3 * V1 - B1*B3;
  FTau = D * (0.5f*VSq + h)*V1 + BSq*V1 - B1*BV;
  FB1 = 0.0f;
  FB2 =  V1*B2 - V2*B1;
  FB3 = -V3*B1 + V1*B3;
  if (DualEnergyFormalism)
    FGE = UGE * V1; 

  /* largest and smallest eigenvectors */
  float ca2 = B1*B1/D;
  float temp1 = cs2 + BSq/D;
  float cf2 = 0.5f * (temp1 + sqrt(fabs(temp1*temp1 - 4.0f*cs2*ca2)));
  float cf = sqrt(cf2);
 
  lp = V1 + cf;
  lm = V1 - cf;
}

// 1D HLL-PLM solver
// We will use the same 1D solver for sweeps in all 3 directions
__forceinline__ __device__ 
void MHD_HLL_PLM_1D(
  float &FD, float &FS1, float &FS2, float &FS3,
  float &FTau, float &FB1, float &FB2,
  float &FB3, float &FPhi, float &FGE,
  float D0, float D1, float D2, float D3,
  float Ein0, float Ein1, float Ein2, float Ein3,
  float V10, float V11, float V12, float V13,
  float V20, float V21, float V22, float V23,
  float V30, float V31, float V32, float V33,
  float B10, float B11, float B12, float B13,
  float B20, float B21, float B22, float B23,
  float B30, float B31, float B32, float B33,
  float Phi0, float Phi1, float Phi2, float Phi3,
  int EOSType, float Ch, float Gamma, float ThetaLimiter, 
  int NeqMHD, int DualEnergyFormalism,
  int i, int j, int k, int idx3d)
{
  float Dl, V1l, V2l, V3l, Einl, B1l, B2l, B3l, Phil;
  float Dr, V1r, V2r, V3r, Einr, B1r, B2r, B3r, Phir;
  float UlD, UlS1, UlS2, UlS3, UlTau, UlB1, UlB2, UlB3, UlPhi, UlGE;
  float UrD, UrS1, UrS2, UrS3, UrTau, UrB1, UrB2, UrB3, UrPhi, UrGE;
  float FlD, FlS1, FlS2, FlS3, FlTau, FlB1, FlB2, FlB3, FlPhi, FlGE;
  float FrD, FrS1, FrS2, FrS3, FrTau, FrB1, FrB2, FrB3, FrPhi, FrGE;
  float lp_l, lp_r, lm_l, lm_r;
  
  PLMPointDevice(Dl  , D0, D1, D2, ThetaLimiter);
  PLMPointDevice(Einl, Ein0, Ein1, Ein2, ThetaLimiter);
  PLMPointDevice(V1l , V10, V11, V12, ThetaLimiter);
  PLMPointDevice(V2l , V20, V21, V22, ThetaLimiter);
  PLMPointDevice(V3l , V30, V31, V32, ThetaLimiter);
  PLMPointDevice(B1l , B10, B11, B12, ThetaLimiter);
  PLMPointDevice(B2l , B20, B21, B22, ThetaLimiter);
  PLMPointDevice(B3l , B30, B31, B32, ThetaLimiter);
  PLMPointDevice(Phil, Phi0, Phi1, Phi2, ThetaLimiter);

  HLLDevice(UlD, UlS1, UlS2, UlS3, UlTau, UlB1, UlB2, UlB3, UlPhi, UlGE,
            FlD, FlS1, FlS2, FlS3, FlTau, FlB1, FlB2, FlB3, FlPhi, FlGE,
            lp_l, lm_l,
            Dl, Einl, V1l, V2l, V3l, B1l, B2l, B3l, Phil, 
            EOSType, Gamma, DualEnergyFormalism, idx3d);

  PLMPointDevice(Dr  , D3, D2, D1, ThetaLimiter);
  PLMPointDevice(Einr, Ein3, Ein2, Ein1, ThetaLimiter);
  PLMPointDevice(V1r , V13, V12, V11, ThetaLimiter);
  PLMPointDevice(V2r , V23, V22, V21, ThetaLimiter);
  PLMPointDevice(V3r , V33, V32, V31, ThetaLimiter);
  PLMPointDevice(B1r , B13, B12, B11, ThetaLimiter);
  PLMPointDevice(B2r , B23, B22, B21, ThetaLimiter);
  PLMPointDevice(B3r , B33, B32, B31, ThetaLimiter);
  PLMPointDevice(Phir, Phi3, Phi2, Phi1, ThetaLimiter);

  HLLDevice(UrD, UrS1, UrS2, UrS3, UrTau, UrB1, UrB2, UrB3, UrPhi, UrGE,
            FrD, FrS1, FrS2, FrS3, FrTau, FrB1, FrB2, FrB3, FrPhi, FrGE,
            lp_r, lm_r,
            Dr, Einr, V1r, V2r, V3r, B1r, B2r, B3r, Phir, 
            EOSType, Gamma, DualEnergyFormalism, idx3d);

  float ap = Max3Device(0.0f, lp_l, lp_r);
  float am = Max3Device(0.0f, -lm_l, -lm_r);

  FD   = (ap*FlD   + am*FrD   - ap*am*(UrD  -UlD  )) / (ap+am);
  FS1  = (ap*FlS1  + am*FrS1  - ap*am*(UrS1 -UlS1 )) / (ap+am);
  FS2  = (ap*FlS2  + am*FrS2  - ap*am*(UrS2 -UlS2 )) / (ap+am);
  FS3  = (ap*FlS3  + am*FrS3  - ap*am*(UrS3 -UlS3 )) / (ap+am);
  FTau = (ap*FlTau + am*FrTau - ap*am*(UrTau-UlTau)) / (ap+am);
  FB1  = (ap*FlB1  + am*FrB1  - ap*am*(UrB1 -UlB1 )) / (ap+am);
  FB2  = (ap*FlB2  + am*FrB2  - ap*am*(UrB2 -UlB2 )) / (ap+am);
  FB3  = (ap*FlB3  + am*FrB3  - ap*am*(UrB3 -UlB3 )) / (ap+am);
  if (DualEnergyFormalism)
    FGE  = (ap*FlGE  + am*FrGE  - ap*am*(UrGE -UlGE )) / (ap+am);
  
  FB1 += UlPhi + 0.5f*(UrPhi-UlPhi) - 0.5f*Ch*(UrB1-UlB1);
  FPhi = UlB1 + 0.5f*(UrB1-UlB1) - 0.5f/Ch*(UrPhi-UlPhi);
  FPhi *= (Ch*Ch);


}

// Kernel for HLL-PLM MHD solver
// 1 thread corresponds to 1 flux at cell interface
__global__ 
void MHD_HLL_PLMKernel(
  float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
  float *FluxB1, float *FluxB2, float *FluxB3, float *FluxPhi, float *FluxGE,
  float *D, float *V1, float *V2, float *V3, float *TE, 
  float *B1, float *B2, float *B3, float *Phi, float *GE,
  int EOSType, float Ch, float ThetaLimiter, float Gamma, int NeqMHD, 
  int DualEnergyFormalism, float min_coeff,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2, int dir)
{
//  int size = dimx*dimy*dimz;
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;

  float D0, D1, D2, D3,
    Ein0, Ein1, Ein2, Ein3,
    V10, V11, V12, V13,
    V20, V21, V22, V23,
    V30, V31, V32, V33,
    B10, B11, B12, B13,
    B20, B21, B22, B23,
    B30, B31, B32, B33,
    Phi0, Phi1, Phi2, Phi3;
  float FD, FS1, FS2, FS3, FTau, FB1, FB2, FB3, FPhi, FGE;

  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);
  
  int iend = i2, jend = j2, kend = k2;
  if (dir == 1) iend = i2+1;
  if (dir == 2) jend = j2+1;
  if (dir == 3) kend = k2+1;
  int idx;
  float Vsq, Bsq;
  if (i >= i1 && i <= iend && j >= j1 && j <= jend && k >= k1 && k <= kend) {
    int m2, m1, p1;
    float *V1R, *V2R, *V3R, *B1R, *B2R, *B3R;

    // Rotate directions to reuse the 1D MHD solver
    if (dir == 1) {
      m2 = -2, m1 = -1, p1 = 1;
      V1R = V1, V2R = V2, V3R = V3;
      B1R = B1, B2R = B2, B3R = B3;
    }
    if (dir == 2) {
      m2 = -2*dimx, m1 = -dimx, p1 = dimx;
      V1R = V2, V2R = V3, V3R = V1;
      B1R = B2, B2R = B3, B3R = B1;
    }
    if (dir == 3) {
      m2 = -2*dimx*dimy, m1 = -dimx*dimy, p1 = dimx*dimy;
      V1R = V3, V2R = V1, V3R = V2;
      B1R = B3, B2R = B1, B3R = B2;
    }

    //
    // Read in the 4 input values for a single flux computation
    //
    idx = idx3d + m2;
    D0 = D[idx];
    V10 = V1R[idx];
    V20 = V2R[idx]; 
    V30 = V3R[idx];
    B10 = B1R[idx];
    B20 = B2R[idx];
    B30 = B3R[idx];
    Phi0 = Phi[idx];
    Vsq = V10*V10 + V20*V20 + V30*V30;
    Bsq = B10*B10 + B20*B20 + B30*B30;
    if (DualEnergyFormalism) 
      Ein0 = GE[idx]; 
    else
      Ein0 = max(TE[idx] - 0.5f*Vsq - 0.5f*Bsq/D0, min_coeff*D0);

    idx = idx3d + m1;
    D1 = D[idx];
    V11 = V1R[idx];
    V21 = V2R[idx]; 
    V31 = V3R[idx];
    B11 = B1R[idx];
    B21 = B2R[idx];
    B31 = B3R[idx];
    Phi1 = Phi[idx];
    Vsq = V11*V11 + V21*V21 + V31*V31;
    Bsq = B11*B11 + B21*B21 + B31*B31;
    if (DualEnergyFormalism)
      Ein1 = GE[idx];
    else 
      Ein1 = max(TE[idx] - 0.5f*Vsq - 0.5f*Bsq/D1, min_coeff*D1);

    idx = idx3d;
    D2 = D[idx];
    V12 = V1R[idx];
    V22 = V2R[idx]; 
    V32 = V3R[idx];
    B12 = B1R[idx];
    B22 = B2R[idx];
    B32 = B3R[idx];
    Phi2 = Phi[idx];
    Vsq = V12*V12 + V22*V22 + V32*V32;
    Bsq = B12*B12 + B22*B22 + B32*B32;
    if (DualEnergyFormalism)
      Ein2 = GE[idx];
    else
      Ein2 = max(TE[idx] - 0.5f*Vsq - 0.5f*Bsq/D2, min_coeff*D2);
    
    idx = idx3d + p1;
    D3 = D[idx];
    V13 = V1R[idx];
    V23 = V2R[idx]; 
    V33 = V3R[idx];
    B13 = B1R[idx];
    B23 = B2R[idx];
    B33 = B3R[idx];
    Phi3 = Phi[idx];
    Vsq = V13*V13 + V23*V23 + V33*V33;
    Bsq = B13*B13 + B23*B23 + B33*B33;
    if (DualEnergyFormalism)
      Ein3 = GE[idx];
    else 
      Ein3 = max(TE[idx] - 0.5f*Vsq - 0.5f*Bsq/D3, min_coeff*D3);

    //
    // Call 1D HLL-PLM solver for flux computation
    //
    MHD_HLL_PLM_1D(
      FD, FS1, FS2, FS3, FTau, 
      FB1, FB2, FB3, FPhi, FGE,
      D0  , D1  , D2  , D3  ,
      Ein0, Ein1, Ein2, Ein3,
      V10 , V11 , V12 , V13 ,
      V20 , V21 , V22 , V23 ,
      V30 , V31 , V32 , V33 ,
      B10 , B11 , B12 , B13 ,
      B20 , B21 , B22 , B23 ,
      B30 , B31 , B32 , B33 ,
      Phi0, Phi1, Phi2, Phi3,
      EOSType, Ch, Gamma, ThetaLimiter, 
      NeqMHD, DualEnergyFormalism, i, j, k, idx3d);

    //
    // Save flux
    //
    FluxD  [idx3d] = FD;
    FluxTau[idx3d] = FTau;
    FluxPhi[idx3d] = FPhi;
    if (DualEnergyFormalism)
      FluxGE[idx3d] = FGE;
    // Rotate directions back
    if (dir == 1) {
      FluxS1[idx3d] = FS1, FluxS2[idx3d] = FS2, FluxS3[idx3d] = FS3;
      FluxB1[idx3d] = FB1, FluxB2[idx3d] = FB2, FluxB3[idx3d] = FB3;
    }
    if (dir == 2) {
      FluxS1[idx3d] = FS3, FluxS2[idx3d] = FS1, FluxS3[idx3d] = FS2;
      FluxB1[idx3d] = FB3, FluxB2[idx3d] = FB1, FluxB3[idx3d] = FB2;
    }
    if (dir == 3) {
      FluxS1[idx3d] = FS2, FluxS2[idx3d] = FS3, FluxS3[idx3d] = FS1;
      FluxB1[idx3d] = FB2, FluxB2[idx3d] = FB3, FluxB3[idx3d] = FB1;
    }

  }
}

extern "C"
void MHD_HLL_PLMGPU(cuMHDData &Data, int dir, float dx)
{
  float min_coeff = 0;
  if (UseMinimumPressureSupport)
    min_coeff = MinimumPressureSupportParameter*0.32*dx*dx/(Gamma*(Gamma-1.0f));

  MHD_HLL_PLMKernel<<<CudaGrid, CudaBlock>>>(
    Data.Flux[IdxD], Data.Flux[IdxS1], Data.Flux[IdxS2], Data.Flux[IdxS3],
    Data.Flux[IdxTau], 
    Data.Flux[IdxB1], Data.Flux[IdxB2], Data.Flux[IdxB3], 
    Data.Flux[IdxPhi], Data.Flux[IdxGE],
    Data.Baryon[IdxD], 
    Data.Baryon[IdxV1], Data.Baryon[IdxV2], Data.Baryon[IdxV3], 
    Data.Baryon[IdxTE], 
    Data.Baryon[IdxB1], Data.Baryon[IdxB2], Data.Baryon[IdxB3], 
    Data.Baryon[IdxPhi], Data.Baryon[IdxGE],
    EOSType, C_h, Theta_Limiter, Gamma, NEQ_MHD, DualEnergyFormalism,
    min_coeff,
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2], dir);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

// Solve for advecting species fields
__global__ void ComputeFluxSpeciesKernel(
  float **FluxSpecies, float **Species, float *FluxD,
  int NSpecies, float ThetaLimiter,
  int dimx, int dimy, int dimz, int i1, int i2, int j1, int j2, int k1, int k2, int dir)
{
//  int size = dimx*dimy*dimz;
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;

  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  int iend = i2, jend = j2, kend = k2;
  if (dir == 1) iend = i2 + 1;
  if (dir == 2) jend = j2 + 1;
  if (dir == 3) kend = k2 + 1;

  if (i >= i1 && i <= iend && j >= j1 && j <= jend && k >= k1 && k <= kend) {
    float s0, s1, s2, s3, s;
    int m2, m1, p1;
    if (dir == 1) m2 = -2, m1 = -1, p1 = 1;
    if (dir == 2) m2 = -2*dimx, m1 = -dimx, p1 = dimx;
    if (dir == 3) m2 = -2*dimx*dimy, m1 = -dimx*dimy, p1 = dimx*dimy;
    for (int f = 0; f < NSpecies; f++) {
      s0 = Species[f][idx3d+m2];
      s1 = Species[f][idx3d+m1];
      s2 = Species[f][idx3d   ];
      s3 = Species[f][idx3d+p1];
    
      if (FluxD[idx3d] >= 0) 
        PLMLeftDevice(s, s0, s1, s2, ThetaLimiter);
      else
        PLMRightDevice(s, s1, s2, s3, ThetaLimiter);
      FluxSpecies[f][idx3d] = s;
    }
    // renormalize
    float norm = 0;
    for (int f = 0; f < NSpecies; f++)
      norm += FluxSpecies[f][idx3d];
    for (int f = 0; f < NSpecies; f++)
      FluxSpecies[f][idx3d] /= norm;
    for (int f = 0; f < NSpecies; f++)
      FluxSpecies[f][idx3d] *= FluxD[idx3d];
  }
}

extern "C"
void ComputeFluxSpeciesGPU(cuMHDData &Data, int dir)
{
  ComputeFluxSpeciesKernel<<<CudaGrid, CudaBlock>>>(
    Data.FluxSpeciesArray, Data.SpeciesArray, Data.Flux[IdxD], 
    NSpecies, Theta_Limiter,
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2], dir);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

// Differentiate fluxes to get dU
__global__ void ComputedUKernel(
  float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau,
  float *dUB1, float *dUB2, float *dUB3, float *dUPhi, float *dUGE,
  float **dUSpecies,
  const float* __restrict FluxD, const float* __restrict FluxS1, 
  const float* __restrict FluxS2, const float* __restrict FluxS3, 
  const float* __restrict FluxTau, const float* __restrict FluxB1, 
  const float* __restrict FluxB2, const float* __restrict FluxB3, 
  const float* __restrict FluxPhi, const float* __restrict FluxGE,
  float **FluxSpecies, int DualEnergyFormalism,
  int NSpecies, float dt, float dx, float Ch,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2, int dir)
{
  const int size = dimx*dimy*dimz;
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);
 
  float dtdx = dt/dx;

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    int iflux = idx3d, ifluxp1;
    if (dir == 1) ifluxp1 = idx3d + 1;
    if (dir == 2) ifluxp1 = idx3d + dimx;
    if (dir == 3) ifluxp1 = idx3d + dimx*dimy;

    dUD  [idx3d] -= (FluxD  [ifluxp1] - FluxD  [iflux]) * dtdx; 
    dUS1 [idx3d] -= (FluxS1 [ifluxp1] - FluxS1 [iflux]) * dtdx;
    dUS2 [idx3d] -= (FluxS2 [ifluxp1] - FluxS2 [iflux]) * dtdx;
    dUS3 [idx3d] -= (FluxS3 [ifluxp1] - FluxS3 [iflux]) * dtdx;
    dUTau[idx3d] -= (FluxTau[ifluxp1] - FluxTau[iflux]) * dtdx;
    dUB1 [idx3d] -= (FluxB1 [ifluxp1] - FluxB1 [iflux]) * dtdx;
    dUB2 [idx3d] -= (FluxB2 [ifluxp1] - FluxB2 [iflux]) * dtdx;
    dUB3 [idx3d] -= (FluxB3 [ifluxp1] - FluxB3 [iflux]) * dtdx;
    dUPhi[idx3d] -= (FluxPhi[ifluxp1] - FluxPhi[iflux]) * dtdx;
    if (DualEnergyFormalism) 
      dUGE [idx3d] -= (FluxGE [ifluxp1] - FluxGE [iflux]) * dtdx;
    for (int f = 0; f < NSpecies; f++)
      dUSpecies[f][idx3d] -= (FluxSpecies[f][ifluxp1] - FluxSpecies[f][iflux]) * dtdx;
  }
}

extern "C"
void ComputedUGPU(cuMHDData &Data, float dt, float dx, int dir)
{
  ComputedUKernel<<<CudaGrid, CudaBlock>>>(
    Data.dU[IdxD], Data.dU[IdxS1], Data.dU[IdxS2], Data.dU[IdxS3],
    Data.dU[IdxTau], 
    Data.dU[IdxB1], Data.dU[IdxB2], Data.dU[IdxB3], 
    Data.dU[IdxPhi], Data.dU[IdxGE],
    Data.dUSpeciesArray,
    Data.Flux[IdxD], Data.Flux[IdxS1], Data.Flux[IdxS2], Data.Flux[IdxS3], 
    Data.Flux[IdxTau], 
    Data.Flux[IdxB1], Data.Flux[IdxB2], Data.Flux[IdxB3], 
    Data.Flux[IdxPhi], Data.Flux[IdxGE],
    Data.FluxSpeciesArray,
    DualEnergyFormalism,
    NSpecies, dt, dx, C_h, 
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2], dir);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void MHDGravitySourceKernel(
  float *dUS1, float *dUS2, float *dUS3, float *dUTau,
  float *D, float *V1, float *V2, float *V3, float dt,
  float *AccelerationField0, 
  float *AccelerationField1, 
  float *AccelerationField2,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    float density = D[idx3d];
    float gx = AccelerationField0[idx3d];
    float gy = AccelerationField1[idx3d];
    float gz = AccelerationField2[idx3d];
    float vx = V1[idx3d];
    float vy = V2[idx3d];
    float vz = V3[idx3d];
    dUS1[idx3d] += dt*gx*density;
    dUS2[idx3d] += dt*gy*density;
    dUS3[idx3d] += dt*gz*density;
    dUTau[idx3d] += dt*density*(gx*vx + gy*vy + gz*vz);
  }
}

extern "C"
void MHDGravitySourceGPU(cuMHDData &Data, float dt)
{
  MHDGravitySourceKernel<<<CudaGrid, CudaBlock>>>(
    Data.dU[IdxS1], Data.dU[IdxS2], Data.dU[IdxS3], Data.dU[IdxTau], 
    Data.Baryon[IdxD], 
    Data.Baryon[IdxV1], Data.Baryon[IdxV2], Data.Baryon[IdxV3], dt,
    Data.AccelerationField[0], 
    Data.AccelerationField[1], 
    Data.AccelerationField[2],
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void MHDDednerSourceKernel(
  float *dUS1, float *dUS2, float *dUS3, float *dUTau,
  float *B1, float *B2, float *B3, float *divB, float *gradPhi,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  const int size = dimx*dimy*dimz;
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    dUS1[idx3d]  -= divB[idx3d] * B1[idx3d];
    dUS2[idx3d]  -= divB[idx3d] * B2[idx3d];
    dUS3[idx3d]  -= divB[idx3d] * B3[idx3d];
    dUTau[idx3d] -= (B1[idx3d]*gradPhi[idx3d] + 
                     B2[idx3d]*gradPhi[size+idx3d] + 
                     B3[idx3d]*gradPhi[2*size+idx3d]);
  }
}
    
    

__global__ 
void MHDDualEnergySourceKernel(
  float *dUGE,
  float *D, float *V1, float *V2, float *V3, float *GE,
  float dtdx, float dtdy, float dtdz, 
  float Gamma, int EOSType, float min_coeff,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    float density = D[idx3d];
    float eint = GE[idx3d];
    eint = max(eint, min_coeff*density);
    int ip1 = idx3d + 1;
    int im1 = idx3d - 1;
    int jp1 = idx3d + dimx;
    int jm1 = idx3d - dimx;
    int kp1 = idx3d + dimx*dimy;
    int km1 = idx3d - dimx*dimy;
    float divVdt = dtdx*(V1[ip1] - V1[im1]) + 
      dtdy*(V2[jp1] - V2[jm1]) + 
      dtdz*(V3[kp1] - V3[km1]);
    float p, h, cs;
    EOSDevice(p, density, eint, h, cs, Gamma, EOSType, 2);
    dUGE[idx3d] -= p*divVdt;
  }
}

extern "C"
void MHDDualEnergySourceGPU(cuMHDData &Data, float dt, float a,
                            float dx, float dy, float dz)
{
  float dtdx = dt/dx/a;
  float dtdy = dt/dy/a;
  float dtdz = dt/dz/a;
  float min_coeff = 0;
  if (UseMinimumPressureSupport)
    min_coeff = MinimumPressureSupportParameter*0.32*dx*dx/(Gamma*(Gamma-1.0f));
  MHDDualEnergySourceKernel<<<CudaGrid, CudaBlock>>>(
    Data.dU[IdxGE],
    Data.Baryon[IdxD], Data.Baryon[IdxV1], Data.Baryon[IdxV2], 
    Data.Baryon[IdxV3], Data.Baryon[IdxGE], 
    dtdx, dtdy, dtdz, Gamma, EOSType, min_coeff,
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}


__global__ 
void MHDComovingSourceKernel(
  float *dUTau, float *dUB1, float *dUB2, float *dUB3, float *dUPhi,
  float *B1, float *B2, float *B3, float dt, float coef,
  int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    float bx = B1[idx3d];
    float by = B2[idx3d];
    float bz = B3[idx3d];
    dUB1[idx3d] += dt*coef*bx;
    dUB2[idx3d] += dt*coef*by;
    dUB3[idx3d] += dt*coef*bz;
    dUTau[idx3d] -= dt*coef*(bx*bx + by*by + bz*bz);
    dUPhi[idx3d] += 0.0f; // Add correct Phi term here ...
  }
}

extern "C"
void MHDComovingSourceGPU(cuMHDData &Data, float dt, float coef)
{
  MHDComovingSourceKernel<<<CudaGrid, CudaBlock>>>(
    Data.dU[IdxTau], Data.dU[IdxB1], Data.dU[IdxB2], Data.dU[IdxB3], 
    Data.dU[IdxPhi], 
    Data.Baryon[IdxB1], Data.Baryon[IdxB2], Data.Baryon[IdxB3], 
    dt, coef,
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void MHDDrivingSourceKernel(
  float *dUS1, float *dUS2, float *dUS3, float *dUTau,
  float *D, float *V1, float *V2, float *V3, 
  float *TE, float *B1, float *B2, float *B3,
  float *DrivingForce1, float *DrivingForce2, 
  float *DrivingForce3, float DrivingEfficiency,
  float dt, int dimx, int dimy, int dimz, 
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    float density = D[idx3d];
    float vx = V1[idx3d];
    float vy = V2[idx3d];
    float vz = V3[idx3d];
    float fx = DrivingForce1[idx3d];
    float fy = DrivingForce2[idx3d];
    float fz = DrivingForce3[idx3d];
    float eint = TE[idx3d] - 0.5f*sqrt(vx*vx + vy*vy + vz*vz) - 
      0.5f*sqrt(pow(B1[idx3d],2)+pow(B2[idx3d],2)+pow(B3[idx3d],2))/density;
    dUS1[idx3d] += dt*density*fx*DrivingEfficiency;
    dUS2[idx3d] += dt*density*fy*DrivingEfficiency;
    dUS3[idx3d] += dt*density*fz*DrivingEfficiency;
    dUTau[idx3d] += dt*density*(fx*vx + fy*vy + fz*vz)*DrivingEfficiency;
  }
}
  
extern "C"
void MHDDrivingSourceGPU(cuMHDData &Data, float dt)
{
  MHDDrivingSourceKernel<<<CudaGrid, CudaBlock>>>(
    Data.dU[IdxS1], Data.dU[IdxS2], Data.dU[IdxS3], Data.dU[IdxTau], 
    Data.Baryon[IdxD], 
    Data.Baryon[IdxV1], Data.Baryon[IdxV2], Data.Baryon[IdxV3], 
    Data.Baryon[IdxTE], 
    Data.Baryon[IdxB1], Data.Baryon[IdxB2], Data.Baryon[IdxB3],
    Data.DrivingForce[0], Data.DrivingForce[1], Data.DrivingForce[2], 
    DrivingEfficiency, dt,
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void Density2FractionKernel(float **Species, float *Density, 
                            float NSpecies, int size)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + blockIdx.x*blockDim.x + 
    threadIdx.x;
  if (idx3d < size) {
    for (int i = 0; i < NSpecies; i++) {
      // float *s = Species[i];
      // s[idx3d] /= Density[idx3d];
      Species[i][idx3d] /= Density[idx3d];
    }
  }
}

extern "C"
void Density2FractionGPU(cuMHDData &Data)
{
  const int size = Data.Dimension[0]*Data.Dimension[1]*Data.Dimension[2];
  Density2FractionKernel<<<CudaGrid, CudaBlock>>>(
    Data.SpeciesArray, Data.Baryon[IdxD], NSpecies, size);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void Fraction2DensityKernel(float **Species, float *Density, 
                            float NSpecies, int size)
{
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + blockIdx.x*blockDim.x + 
    threadIdx.x;
  if (idx3d < size) {
    for (int i = 0; i < NSpecies; i++) {
      Species[i][idx3d] *= Density[idx3d];
    }
  }
}

extern "C"
void Fraction2DensityGPU(cuMHDData &Data)
{
  const int size = Data.Dimension[0]*Data.Dimension[1]*Data.Dimension[2];
  Fraction2DensityKernel<<<CudaGrid, CudaBlock>>>(
    Data.SpeciesArray, Data.Baryon[IdxD], NSpecies, size);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void UpdateMHDPrimKernel(
  float *D, float *V1, float *V2, float *V3, float *TE,
  float *B1, float *B2, float *B3, float *Phi, float *GE, float **Species,
  float *OldD, float *OldV1, float *OldV2, float *OldV3, float *OldTE,
  float *OldB1, float *OldB2, float *OldB3, float *OldPhi, 
  float *OldGE, float **OldSpecies,
  float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau,
  float *dUB1, float *dUB2, float *dUB3, float *dUPhi, 
  float *dUGE, float **dUSpecies,
  int DualEnergyFormalism, int DualEnergyFormalismEta1,
  float Gamma, int EOSType, float SmallT, float Mu,
  int NSpecies, int RKStep, float dt, float Ch, float Cp,
  int dimx, int dimy, int dimz,
  int i1, int i2, int j1, int j2, int k1, int k2)
{
  //const int size = dimx*dimy*dimz;
  float D_n, V1_n, V2_n, V3_n, TE_n, B1_n, B2_n, B3_n, Phi_n, GE_n,
    Tau_n, S1_n, S2_n, S3_n;
  float D_1, V1_1, V2_1, V3_1, TE_1, B1_1, B2_1, B3_1, Phi_1, Tau_1, GE_1,
    S1_1, S2_1, S3_1;
  float D_np1, B1_np1, B2_np1, B3_np1, Phi_np1, Tau_np1, GE_np1,
    S1_np1, S2_np1, S3_np1;
  int idx3d = blockIdx.y*blockDim.x*gridDim.x + 
    blockIdx.x*blockDim.x + threadIdx.x + i1 + j1*dimx + k1*dimx*dimy;
  
  int i = idx3d % dimx;
  int j = (idx3d % (dimx*dimy)) / dimx;
  int k = idx3d / (dimx*dimy);

  if (i >= i1 && i <= i2 && j >= j1 && j <= j2 && k >= k1 && k <= k2) {
    // update species
    float sp, sum = 0, norm = 0;
    for (int f = 0; f < NSpecies; f++) {
      if (RKStep == 1) {
        sp = Species[f][idx3d]*D[idx3d] + dUSpecies[f][idx3d];
      } else {
        sp = 0.5f*(OldSpecies[f][idx3d] + Species[f][idx3d]*D[idx3d] + 
                   dUSpecies[f][idx3d]);
      }
      Species[f][idx3d] = sp;
      sum += sp;
    }
    for (int f = 0; f < NSpecies; f++) {
      Species[f][idx3d] = min(1.0f, max(Species[f][idx3d]/sum, 1e-20));
      Species[f][idx3d] = Species[f][idx3d] / sum;
      norm += Species[f][idx3d];
    }
    for (int f = 0; f < NSpecies; f++) 
      Species[f][idx3d] /= norm;

    D_n   = OldD[idx3d];
    V1_n  = OldV1[idx3d];
    V2_n  = OldV2[idx3d];
    V3_n  = OldV3[idx3d];
    TE_n  = OldTE[idx3d];
    B1_n  = OldB1[idx3d];
    B2_n  = OldB2[idx3d];
    B3_n  = OldB3[idx3d];
    Phi_n = OldPhi[idx3d];
    if (DualEnergyFormalism)
      GE_n = OldGE[idx3d];
      
    Tau_n = D_n * TE_n;
    S1_n  = D_n * V1_n;
    S2_n  = D_n * V2_n;
    S3_n  = D_n * V3_n;
    if (RKStep == 1) {
      D_np1   = D_n   + dUD  [idx3d];
      S1_np1  = S1_n  + dUS1 [idx3d];
      S2_np1  = S2_n  + dUS2 [idx3d];
      S3_np1  = S3_n  + dUS3 [idx3d];
      Tau_np1 = Tau_n + dUTau[idx3d];
      B1_np1  = B1_n  + dUB1 [idx3d];
      B2_np1  = B2_n  + dUB2 [idx3d];
      B3_np1  = B3_n  + dUB3 [idx3d];
      Phi_np1 = Phi_n + dUPhi[idx3d];
      if (DualEnergyFormalism) 
        GE_np1  = D_n*GE_n + dUGE[idx3d];
    } else {
      D_1   = D[idx3d];
      V1_1  = V1[idx3d];
      V2_1  = V2[idx3d];
      V3_1  = V3[idx3d];
      TE_1  = TE[idx3d];
      B1_1  = B1[idx3d];
      B2_1  = B2[idx3d];
      B3_1  = B3[idx3d];
      Phi_1 = Phi[idx3d];
      if (DualEnergyFormalism)
        GE_1  = GE[idx3d];

      Tau_1 = D_1 * TE_1;
      S1_1  = D_1 * V1_1;
      S2_1  = D_1 * V2_1;
      S3_1  = D_1 * V3_1;

      D_np1   = 0.5f*(D_n   + D_1   + dUD  [idx3d]);
      S1_np1  = 0.5f*(S1_n  + S1_1  + dUS1 [idx3d]);
      S2_np1  = 0.5f*(S2_n  + S2_1  + dUS2 [idx3d]);
      S3_np1  = 0.5f*(S3_n  + S3_1  + dUS3 [idx3d]);
      Tau_np1 = 0.5f*(Tau_n + Tau_1 + dUTau[idx3d]);
      B1_np1  = 0.5f*(B1_n  + B1_1  + dUB1 [idx3d]);
      B2_np1  = 0.5f*(B2_n  + B2_1  + dUB2 [idx3d]);
      B3_np1  = 0.5f*(B3_n  + B3_1  + dUB3 [idx3d]);
      Phi_np1 = 0.5f*(Phi_n + Phi_1 + dUPhi[idx3d]);
      if (DualEnergyFormalism)
        GE_np1  = 0.5f*(D_n*GE_n + D_1*GE_1 + dUGE[idx3d]);
    }

    if (D_np1 < 0 || isnan(D_np1))
      printf("UpdateMHDPrimKernel: rho < 0 at (%d,%d,%d)\n", i, j, k);

    float vx = S1_np1 / D_np1;
    float vy = S2_np1 / D_np1;
    float vz = S3_np1 / D_np1;
    float etot = Tau_np1 / D_np1;

    D  [idx3d] = D_np1;
    V1 [idx3d] = vx;
    V2 [idx3d] = vy;
    V3 [idx3d] = vz;
    TE [idx3d] = etot;
    B1 [idx3d] = B1_np1;
    B2 [idx3d] = B2_np1;
    B3 [idx3d] = B3_np1;
    Phi[idx3d] = Phi_np1 * exp(-dt*(Ch/Cp)*(Ch/Cp));

    if (DualEnergyFormalism) {
      float v2 = vx*vx + vy*vy + vz*vz;
      float B2 = B1_np1*B1_np1 + B2_np1*B2_np1 + B3_np1*B3_np1;

      float eint = GE_np1 / D_np1;
      float eint1 = etot - 0.5f*v2 - 0.5f*B2/D_np1;
      float p, h, cs;
      if (eint1 > 0) 
        EOSDevice(p, D_np1, eint1, h, cs, Gamma, EOSType, 2);
      else
        cs = 0.0;
      if (cs*cs > DualEnergyFormalismEta1*v2 &&
          cs*cs > DualEnergyFormalismEta1*B2/D_np1 &&
          eint1 > 0.5f*eint)
        eint = eint1;
      float emin = SmallT/(Mu*(Gamma-1.0f));
      eint = max(eint, emin);
      GE[idx3d] = eint;
      TE[idx3d] = eint + 0.5f*v2 + 0.5f*B2/D_np1;
    }
  }
}

extern "C"
void UpdateMHDPrimGPU(cuMHDData &Data, int RKStep, float dt)
{
  UpdateMHDPrimKernel<<<CudaGrid, CudaBlock>>>(
    Data.Baryon[IdxD], 
    Data.Baryon[IdxV1], Data.Baryon[IdxV2], Data.Baryon[IdxV3], 
    Data.Baryon[IdxTE], 
    Data.Baryon[IdxB1], Data.Baryon[IdxB2], Data.Baryon[IdxB3],
    Data.Baryon[IdxPhi], Data.Baryon[IdxGE], Data.SpeciesArray,
    Data.OldBaryon[IdxD], 
    Data.OldBaryon[IdxV1], Data.OldBaryon[IdxV2], Data.OldBaryon[IdxV3],
    Data.OldBaryon[IdxTE], 
    Data.OldBaryon[IdxB1], Data.OldBaryon[IdxB2], Data.OldBaryon[IdxB3],
    Data.OldBaryon[IdxPhi], Data.OldBaryon[IdxGE], Data.OldSpeciesArray,
    Data.dU[IdxD], 
    Data.dU[IdxS1], Data.dU[IdxS2], Data.dU[IdxS3], Data.dU[IdxTau], 
    Data.dU[IdxB1], Data.dU[IdxB2], Data.dU[IdxB3], Data.dU[IdxPhi], 
    Data.dU[IdxGE], Data.dUSpeciesArray,
    DualEnergyFormalism, DualEnergyFormalismEta1,
    Gamma, EOSType, SmallT, Mu,
    NSpecies, RKStep, dt, C_h, C_p, 
    Data.Dimension[0], Data.Dimension[1], Data.Dimension[2],
    Data.StartIndex[0], Data.EndIndex[0],
    Data.StartIndex[1], Data.EndIndex[1],
    Data.StartIndex[2], Data.EndIndex[2]);
  CUDA_SAFE_CALL( cudaGetLastError() );
}

__global__ 
void MHDSaveSubgridFluxKernel(
  float *LeftFluxes,  float *RightFluxes, const float* __restrict Flux,
  float dtdx, int dimx, int dimy, int dimz,
  int fistart, int fiend, int fjstart, int fjend, 
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
    if (dir == 1) {
      lindex = i*dimx + j*dimx*dimy + lface;
      rindex = i*dimx + j*dimx*dimy + rface;
    } else if (dir == 2) {
      lindex = i + j*dimx*dimy + lface*dimx;
      rindex = i + j*dimx*dimy + rface*dimx;
    } else if (dir == 3) {
      lindex = i + j*dimx + lface*dimx*dimy;
      rindex = i + j*dimx + rface*dimx*dimy; 
    }
    LeftFluxes[offset] = 0.5f*Flux[lindex]*dtdx;
    RightFluxes[offset] = 0.5f*Flux[rindex]*dtdx;
  }
}

extern "C"
void MHDSaveSubgridFluxGPU(float *LeftFlux, float *RightFlux, float *Flux,
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
  CUDA_SAFE_CALL( cudaGetLastError() );
  for (int i = 0; i < nf; i++) {
    LeftFlux[i] += Tmp3[i];
    RightFlux[i] += Tmp4[i];
  }
}                     
