#ifndef _CUMHD_H_
#define _CUMHD_H_

#define NMHD 10
#define MAX_SPECIES 20
#define NBLOCK 128
#define PLM_GHOST_SIZE 2

enum {IdxD, IdxV1, IdxV2, IdxV3, IdxTE, IdxB1, IdxB2, IdxB3, IdxPhi, IdxGE};
enum {IdxS1 = 1, IdxS2 = 2, IdxS3 = 3, IdxTau = 4};

typedef struct {
  // baryon fields (primitives)
  float *Baryon[NMHD];
  // old baryon fields
  float *OldBaryon[NMHD];
  // fluxes
  float *Flux[NMHD];
  // dU
  float *dU[NMHD];
  // source terms
  float  *AccelerationField[3], *DrivingForce[3];
  // Color
  float *Species[MAX_SPECIES];
  float *OldSpecies[MAX_SPECIES];
  float *FluxSpecies[MAX_SPECIES];
  float *dUSpecies[MAX_SPECIES];
  float **SpeciesArray;
  float **OldSpeciesArray;
  float **FluxSpeciesArray;
  float **dUSpeciesArray;
  // Parameters
  int Dimension[3], StartIndex[3], EndIndex[3];
  void *GPUMem;
  size_t GPUMemSize;
  size_t GPUMemOffset;
} cuMHDData;

extern "C"
void MHD_HLL_PLMGPU(cuMHDData &Data,  int dir, float dx);

extern "C"
void ComputeFluxSpeciesGPU(cuMHDData &Data, int dir);

extern "C"
void ComputedUGPU(cuMHDData &Data, float dt, float dx, int dir);

extern "C"
void MHDGravitySourceGPU(cuMHDData &Data, float dt);

extern "C"
void MHDComovingSourceGPU(cuMHDData &Data, float dt, float coef);

extern "C"
void MHDDrivingSourceGPU(cuMHDData &Data, float dt);

extern "C"
void MHDDualEnergySourceGPU(cuMHDData &Data, float dt, float a,
                            float dx, float dy, float dz);

extern "C"
void Density2FractionGPU(cuMHDData &Data);

extern "C"
void Fraction2DensityGPU(cuMHDData &Data);

extern "C"
void UpdateMHDPrimGPU(cuMHDData &Data, int RKStep, float dt);

extern "C"
void MHDSaveSubgridFluxGPU(float *LeftFlux, float *RightFlux, float *Flux,
                           float *Tmp1, float *Tmp2, float *Tmp3, float *Tmp4, 
                           float dtdx,
                           const int dimx, const int dimy, const int dimz,
                           const int fistart, const int fiend,
                           const int fjstart, const int fjend,
                           const int lface, const int rface, int dir);
#endif
