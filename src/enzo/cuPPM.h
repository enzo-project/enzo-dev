#ifndef _CUPPM_H_
#define _CUPPM_H_

#include "CUDAUtil.h"
#include "fortran.def"
#define tiny 1e-20

typedef struct {
  // Input variables
  float *d, *e, *p, *u, *v, *w, *ge, *col, *gr[3];
  // Transposed storage
  float *dslice, *eslice, *pslice, *uslice, *vslice, *wslice, 
    *geslice, *colslice, *grslice;
  // Flux variables
  float *df, *ef, *uf, *vf, *wf, *gef, *colf;
  // PPM temporaray variables
  float *dls, *drs, *pls, *prs, *gels, *gers, *uls, *urs,
    *vls, *vrs, *wls, *wrs, *pbar, *ubar, *colls, *colrs, *ges, *gesf,
    *diffcoef, *flatten;
  int *fallback;
} cuPPMData;

typedef struct {
  // BaryonField index
  int DensNum, TENum, Vel1Num, Vel2Num, Vel3Num, GENum;
  int MaxNumberOfBaryonFields;
  // Grid parameters
  Elong_int GridGlobalStart[3];
  int GridStartIndex[3];
  int GridEndIndex[3];
  int GridDimension[3];
  float dx[3];
  // Hydro parameters
  float Gamma;
  int PPMDiffusionParameter;
  int PPMFlatteningParameter;
  int PPMSteepeningParameter;
  int ConservativeReconstruction;
  int PositiveReconstruction;
  int GravityOn;
  int DualEnergyFormalism;
  float DualEnergyFormalismEta1;
  float DualEnergyFormalismEta2;
  int PressureFree;
  float MinimumPressure;
  int RiemannSolver;
  int RiemannSolverFallback;
  // Color parameters
  int NumberOfColours;
  int colnum[MAX_COLOR];
} cuPPMParameter;

extern "C" 
void cuPPMInitParameter(cuPPMParameter *Para,
                        int DensNum, int TENum,
                        int Vel1Num, int Vel2Num, int Vel3Num,
                        int GENum, int MaxNumberOfBaryonFields,
                        Elong_int GridGlobalStart[],
                        int GridStartIndex[],
                        int GridEndIndex[],
                        int GridDimension[], 
                        float dxx, float dxy, float dxz,
                        float Gamma,
                        int PPMDiffusionParameter,
                        int PPMFlatteningParameter,
                        int PPMSteepeningParameter,
                        int ConservativeReconstruction,
                        int PositiveReconstruction,
                        int GravityOn,
                        int DualEnergyFormalism,
                        float DualEnergyFormalismEta1,
                        float DualEnergyFormalismEta2,
                        int PressureFree,
                        float MinimumPressure,
                        int RiemannSolver,
                        int RiemannSolverFallback,
                        int NumberOfColours, int *colnum);

extern "C"
void cuPPMInitData(cuPPMData *Data,
                   cuPPMParameter &Para);

extern "C"
void cuPPMDestroy(cuPPMData &Data, cuPPMParameter &Para);

extern "C"
void cuPPMSetBaryon(cuPPMData &Data, cuPPMParameter &Para,
                    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS],
                    float *Pressure,
                    float *AccelerationField[MAX_DIMENSION]);

extern "C"
void cuPPMGetBaryon(cuPPMData &Data, cuPPMParameter &Para,
                    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS]);

extern "C"
void cuPPMSweep(cuPPMData &Data, cuPPMParameter &Para, float dtFixed, int dir);

extern "C"
void cuPPMSaveSubgridFluxes(
  fluxes *SubgridFluxes[], int NumberOfSubgrids, Elong_int GridGlobalStart[],
  cuPPMData &Data, cuPPMParameter &Para, int dir);

#endif
