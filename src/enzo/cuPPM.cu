/***********************************************************************
/
/  SOLVE HYDRO EQUATIONS USING PPM METHOD ON GPU
/
/  written by: Peng Wang 
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/         Solve hydro equations using PPM method on GPU
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "Fluxes.h"
#include "CUDAUtil.h"
#include "SaveSubgridFluxCUDA.h"
#include "cuPPM.h"
#include "cuPPM.cuh" 

//
// Initializes various parameters the cuPPM solver needs
//
extern "C" 
void cuPPMInitParameter(cuPPMParameter *Para,
                        int DensNum, int TENum,
                        int Vel1Num, int Vel2Num, int Vel3Num, int GENum,
                        int MaxNumberOfBaryonFields,
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
                        int NumberOfColours, int *colnum)
{
  if (GridDimension[0] <= 1 || 
      GridDimension[1] <= 1 ||
      GridDimension[2] <= 1) {
    printf("GPU PPM solver works only for 3D.\n");
    exit(1);
  }
  if (PPMDiffusionParameter) {
    printf("PPM Diffusion not supported on GPU.\n");
    exit(1);
  }
  if (ConservativeReconstruction) {
    printf("ConvervativeRecontruction not supported on GPU.\n");
    exit(1);
  }
  if (PositiveReconstruction) {
    printf("PositiveReconstruction not supported on GPU.\n");
    exit(1);
  }
  Para->DensNum = DensNum;
  Para->TENum = TENum;
  Para->Vel1Num = Vel1Num;
  Para->Vel2Num = Vel2Num;
  Para->Vel3Num = Vel3Num;
  Para->GENum = GENum;
  Para->MaxNumberOfBaryonFields = MaxNumberOfBaryonFields;
  for (int i = 0; i < 3; i++) {
    Para->GridGlobalStart[i] = GridGlobalStart[i];
    Para->GridStartIndex[i] = GridStartIndex[i];
    Para->GridEndIndex[i] = GridEndIndex[i];
    Para->GridDimension[i] = GridDimension[i];
  }
  Para->dx[0] = dxx;
  Para->dx[1] = dxy;
  Para->dx[2] = dxz;
  Para->Gamma = Gamma;
  Para->PPMDiffusionParameter = PPMDiffusionParameter;
  Para->PPMFlatteningParameter = PPMFlatteningParameter;
  Para->PPMSteepeningParameter = PPMSteepeningParameter;
  Para->ConservativeReconstruction = ConservativeReconstruction;
  Para->PositiveReconstruction = PositiveReconstruction;
  Para->GravityOn = GravityOn;
  Para->DualEnergyFormalism = DualEnergyFormalism;
  Para->DualEnergyFormalismEta1 = DualEnergyFormalismEta1;
  Para->DualEnergyFormalismEta2 = DualEnergyFormalismEta2;
  Para->PressureFree = PressureFree;
  Para->MinimumPressure = MinimumPressure;
  Para->RiemannSolver = RiemannSolver;
  Para->RiemannSolverFallback = RiemannSolverFallback;
  Para->NumberOfColours = NumberOfColours;
  for (int i = 0; i < NumberOfColours; i++)
    Para->colnum[i] = colnum[i];
}
//
// Allocate arrays on GPU given the parameters
//                              
extern "C"
void cuPPMInitData(cuPPMData *Data, cuPPMParameter &Para)
{
  size_t size = 1;
  for (int i = 0; i < 3; i++) {
    size *= Para.GridDimension[i];
  }
  const size_t sizebytes = size*sizeof(float);
  // input variables
  cumalloc((void**)&Data->d, sizebytes);
  cumalloc((void**)&Data->e, sizebytes);
  cumalloc((void**)&Data->p, sizebytes);
  cumalloc((void**)&Data->u, sizebytes);
  cumalloc((void**)&Data->v, sizebytes);
  cumalloc((void**)&Data->w, sizebytes);
  if (Para.GravityOn) 
    for (int i = 0; i < 3; i++)
      cumalloc((void**)&Data->gr[i], sizebytes);
  if (Para.DualEnergyFormalism) 
    cumalloc((void**)&Data->ge, sizebytes);
  if (Para.NumberOfColours > 0) 
    cumalloc((void**)&Data->col, Para.NumberOfColours*sizebytes);
  // transposed storage
  cumalloc((void**)&Data->dslice, sizebytes);
  cumalloc((void**)&Data->eslice, sizebytes);
  cumalloc((void**)&Data->pslice, sizebytes);
  cumalloc((void**)&Data->uslice, sizebytes);
  cumalloc((void**)&Data->vslice, sizebytes);
  cumalloc((void**)&Data->wslice, sizebytes);
  if (Para.GravityOn) 
    cumalloc((void**)&Data->grslice, sizebytes);
  if (Para.DualEnergyFormalism) 
    cumalloc((void**)&Data->geslice, sizebytes);
  if (Para.NumberOfColours > 0) 
    cumalloc((void**)&Data->colslice, Para.NumberOfColours*sizebytes);
  // fluxes
  cumalloc((void**)&Data->df, sizebytes);
  cumalloc((void**)&Data->ef, sizebytes);
  cumalloc((void**)&Data->uf, sizebytes);
  cumalloc((void**)&Data->vf, sizebytes);
  cumalloc((void**)&Data->wf, sizebytes);
  if (Para.DualEnergyFormalism) 
    cumalloc((void**)&Data->gef, sizebytes);
  if (Para.NumberOfColours > 0) 
    cumalloc((void**)&Data->colf, Para.NumberOfColours*sizebytes);
  // temporaries
  cumalloc((void**)&Data->dls, sizebytes);
  cumalloc((void**)&Data->drs, sizebytes);
  cumalloc((void**)&Data->pls, sizebytes);
  cumalloc((void**)&Data->prs, sizebytes);
  cumalloc((void**)&Data->gels, sizebytes);
  cumalloc((void**)&Data->gers, sizebytes);
  cumalloc((void**)&Data->uls, sizebytes);
  cumalloc((void**)&Data->urs, sizebytes);
  cumalloc((void**)&Data->vls, sizebytes);
  cumalloc((void**)&Data->vrs, sizebytes);
  cumalloc((void**)&Data->wls, sizebytes);
  cumalloc((void**)&Data->wrs, sizebytes);
  cumalloc((void**)&Data->pbar, sizebytes);
  cumalloc((void**)&Data->ubar, sizebytes);
  if (Para.DualEnergyFormalism) {
    cumalloc((void**)&Data->ges, sizebytes);
    cumalloc((void**)&Data->gesf, sizebytes);
  }
  if (Para.NumberOfColours > 0) {
    cumalloc((void**)&Data->colls, Para.NumberOfColours*sizebytes);
    cumalloc((void**)&Data->colrs, Para.NumberOfColours*sizebytes);
  } 
  if (Para.PPMDiffusionParameter) 
    cumalloc((void**)&Data->diffcoef, sizebytes);
  if (Para.PPMFlatteningParameter) 
    cumalloc((void**)&Data->flatten, sizebytes);
  if (Para.RiemannSolverFallback) {
    cumalloc((void**)&Data->fallback, sizeof(int));
    cudaMemset(Data->fallback, 0, sizeof(int));
  }
}

//
// Free space on GPU
//
extern "C"
void cuPPMDestroy(cuPPMData &Data, cuPPMParameter &Para)
{
  // input variables
  cufree(Data.d);
  cufree(Data.e);
  cufree(Data.p);
  cufree(Data.u);
  cufree(Data.v);
  cufree(Data.w);
  if (Para.GravityOn) 
    for (int i = 0; i < 3; i++)
      cufree(Data.gr[i]);
  if (Para.DualEnergyFormalism)
    cufree(Data.ge);
  if (Para.NumberOfColours > 0)
    cufree(Data.col);
  // transposed storage
  cufree(Data.dslice);
  cufree(Data.eslice);
  cufree(Data.pslice);
  cufree(Data.uslice);
  cufree(Data.vslice);
  cufree(Data.wslice);
  if (Para.GravityOn)
    cufree(Data.grslice);
  if (Para.DualEnergyFormalism)
    cufree(Data.geslice);
  if (Para.NumberOfColours > 0)
    cufree(Data.colslice);
  // fluxes
  cufree(Data.df);
  cufree(Data.ef);
  cufree(Data.uf);
  cufree(Data.vf);
  cufree(Data.wf);
  if (Para.DualEnergyFormalism) 
    cufree(Data.gef);
  if (Para.NumberOfColours > 0)
    cufree(Data.colf);
  // temporaries
  cufree(Data.dls);
  cufree(Data.drs);
  cufree(Data.pls);
  cufree(Data.prs);
  cufree(Data.gels);
  cufree(Data.gers);
  cufree(Data.uls);
  cufree(Data.urs);
  cufree(Data.vls);
  cufree(Data.vrs);
  cufree(Data.wls);
  cufree(Data.wrs);
  cufree(Data.pbar);
  cufree(Data.ubar);
  if (Para.DualEnergyFormalism) {
    cufree(Data.ges);
    cufree(Data.gesf);
  }
  if (Para.NumberOfColours > 0) {
    cufree(Data.colls);
    cufree(Data.colrs);
  }
  if (Para.PPMDiffusionParameter)
    cufree(Data.diffcoef);
  if (Para.PPMFlatteningParameter)
    cufree(Data.flatten);
  if (Para.RiemannSolverFallback)
    cufree(Data.fallback);
}
    
//
// Transfer all BaryonFields from CPU to GPU
//            
extern "C"
void cuPPMSetBaryon(cuPPMData &Data, cuPPMParameter &Para,
                    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS],
                    float *Pressure,
                    float *AccelerationField[MAX_DIMENSION])
{
  const size_t size = 
    Para.GridDimension[0]*Para.GridDimension[1]*Para.GridDimension[2];
  const size_t sizebytes = size*sizeof(float);
  cudaMemcpy(Data.d, BaryonField[Para.DensNum], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(Data.e, BaryonField[Para.TENum], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(Data.p, Pressure, sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(Data.u, BaryonField[Para.Vel1Num], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(Data.v, BaryonField[Para.Vel2Num], sizebytes, cudaMemcpyHostToDevice);
  cudaMemcpy(Data.w, BaryonField[Para.Vel3Num], sizebytes, cudaMemcpyHostToDevice);
  if (Para.GravityOn) {
    for (int i = 0; i < 3; i++) 
      cudaMemcpy(Data.gr[i], AccelerationField[i], sizebytes, 
                 cudaMemcpyHostToDevice);
  }
  if (Para.DualEnergyFormalism)
    cudaMemcpy(Data.ge, BaryonField[Para.GENum], sizebytes, cudaMemcpyHostToDevice);
  if (Para.NumberOfColours > 0) {
    for (int c = 0; c < Para.NumberOfColours; c++)
      cudaMemcpy(Data.col+c*size, BaryonField[Para.colnum[c]], sizebytes,
                 cudaMemcpyHostToDevice);
  }
}
//
// Transfer all BaryonFields from GPU to CPU
//
extern "C"
void cuPPMGetBaryon(cuPPMData &Data, cuPPMParameter &Para,
                    float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS])
{
  const size_t size = 
    Para.GridDimension[0]*Para.GridDimension[1]*Para.GridDimension[2];
  const size_t sizebytes = size*sizeof(float);
  cudaMemcpy(BaryonField[Para.DensNum], Data.d, sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Para.TENum], Data.e, sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Para.Vel1Num], Data.u, sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Para.Vel2Num], Data.v, sizebytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(BaryonField[Para.Vel3Num], Data.w, sizebytes, cudaMemcpyDeviceToHost);
  if (Para.DualEnergyFormalism)
    cudaMemcpy(BaryonField[Para.GENum], Data.ge, sizebytes, cudaMemcpyDeviceToHost);
  if (Para.NumberOfColours > 0) {
    for (int i = 0; i < Para.NumberOfColours; i++)
      cudaMemcpy(BaryonField[Para.colnum[i]], Data.col+i*size, size*sizeof(float),
                 cudaMemcpyDeviceToHost);
  }
}

//
// PPM sweep on GPU
//
extern "C"
void cuPPMSweep(cuPPMData &Data, cuPPMParameter &Para, float dtFixed, int dir)
{
  //
  // Prepare parameter for calling 1D PPM solver
  //
  const int dimx = Para.GridDimension[0];
  const int dimy = Para.GridDimension[1];
  const int dimz = Para.GridDimension[2];
  const int size = dimx * dimy * dimz;
  int idim, jdim, kdim, i1, i2, j1, j2, k1, k2;
  if (dir == 0) {
    i1 = Para.GridStartIndex[0], i2 = Para.GridEndIndex[0];
    j1 = 0, j2 = Para.GridDimension[1] - 1;
    k1 = 0, k2 = Para.GridDimension[2] - 1;
    idim = Para.GridDimension[0];
    jdim = Para.GridDimension[1];
    kdim = Para.GridDimension[2];
  } else if (dir == 1) {
    i1 = Para.GridStartIndex[1], i2 = Para.GridEndIndex[1];
    j1 = 0, j2 = Para.GridDimension[2] - 1;
    k1 = 0, k2 = Para.GridDimension[0] - 1;
    idim = Para.GridDimension[1];
    jdim = Para.GridDimension[2];
    kdim = Para.GridDimension[0];
  } else if (dir == 2) {
    i1 = Para.GridStartIndex[2], i2 = Para.GridEndIndex[2];
    j1 = 0, j2 = Para.GridDimension[0] - 1;
    k1 = 0, k2 = Para.GridDimension[1] - 1;
    idim = Para.GridDimension[2];
    jdim = Para.GridDimension[0];
    kdim = Para.GridDimension[1];
  }
  float dx = Para.dx[dir];

  //
  // Copy data into transpoed storage
  //  This trick allows using the same 1D PPM solver for all 3 directions
  //
  cuPPMData Data1;
  memcpy(&Data1, &Data, sizeof(cuPPMData));
  if (dir == 0) {
    // in 1D slice data is the same as baryonfield, so just pointer copy
    Data1.dslice = Data.d, Data1.pslice = Data.p, Data1.eslice = Data.e;
    Data1.uslice = Data.u, Data1.vslice = Data.v, Data1.wslice = Data.w;
    Data1.grslice = Data.gr[0], Data1.geslice = Data.ge;
    Data1.colslice = Data.col;
  } else if (dir == 1) {
    Transpose<1>(Data1.dslice, Data.d, dimx, dimy, dimz);
    Transpose<1>(Data1.eslice, Data.e, dimx, dimy, dimz);
    Transpose<1>(Data1.pslice, Data.p, dimx, dimy, dimz);
    Transpose<1>(Data1.uslice, Data.v, dimx, dimy, dimz);
    Transpose<1>(Data1.vslice, Data.w, dimx, dimy, dimz);
    Transpose<1>(Data1.wslice, Data.u, dimx, dimy, dimz);
    if (Para.GravityOn)
      Transpose<1>(Data1.grslice, Data.gr[1], dimx, dimy, dimz);
    if (Para.DualEnergyFormalism)
      Transpose<1>(Data1.geslice, Data.ge, dimx, dimy, dimz);
    if (Para.NumberOfColours > 0)
      for (int c = 0; c < Para.NumberOfColours; c++)
        Transpose<1>(Data1.colslice+c*size, Data.col+c*size, dimx, dimy, dimz);
  } else if (dir == 2) {
    Transpose<2>(Data1.dslice, Data.d, dimx, dimy, dimz);
    Transpose<2>(Data1.eslice, Data.e, dimx, dimy, dimz);
    Transpose<2>(Data1.pslice, Data.p, dimx, dimy, dimz);
    Transpose<2>(Data1.uslice, Data.w, dimx, dimy, dimz);
    Transpose<2>(Data1.vslice, Data.u, dimx, dimy, dimz);
    Transpose<2>(Data1.wslice, Data.v, dimx, dimy, dimz);
    if (Para.GravityOn)
      Transpose<2>(Data1.grslice, Data.gr[2], dimx, dimy, dimz);
    if (Para.DualEnergyFormalism)
      Transpose<2>(Data1.geslice, Data.ge, dimx, dimy, dimz);
    if (Para.NumberOfColours > 0)
      for (int c = 0; c < Para.NumberOfColours; c++)
        Transpose<2>(Data1.colslice+c*size, Data.col+c*size, dimx, dimy, dimz);
  } 
  //
  // Call PPM solver on CUDA
  //
  PPM_CUDA(Data1, Para, idim, jdim, kdim, i1, i2, j1, j2, k1, k2, dx, dtFixed);
  //
  // Tranpose data layout back
  //
  if (dir == 1) {
    Transpose<2>(Data.d, Data1.dslice, idim, jdim, kdim);
    Transpose<2>(Data.e, Data1.eslice, idim, jdim, kdim);
    Transpose<2>(Data.v, Data1.uslice, idim, jdim, kdim);
    Transpose<2>(Data.w, Data1.vslice, idim, jdim, kdim);
    Transpose<2>(Data.u, Data1.wslice, idim, jdim, kdim);
    if (Para.DualEnergyFormalism)
      Transpose<2>(Data.ge, Data.geslice, idim, jdim, kdim);
    if (Para.NumberOfColours > 0)
      for (int c = 0; c < Para.NumberOfColours; c++)
        Transpose<2>(Data.col+c*size, Data.colslice+c*size, idim, jdim, kdim);
  } else if (dir == 2) {
    Transpose<1>(Data.d, Data1.dslice, idim, jdim, kdim);
    Transpose<1>(Data.e, Data1.eslice, idim, jdim, kdim);
    Transpose<1>(Data.w, Data1.uslice, idim, jdim, kdim);
    Transpose<1>(Data.u, Data1.vslice, idim, jdim, kdim);
    Transpose<1>(Data.v, Data1.wslice, idim, jdim, kdim);
    if (Para.DualEnergyFormalism)
      Transpose<1>(Data.ge, Data.geslice, idim, jdim, kdim);
    if (Para.NumberOfColours > 0)
      for (int c = 0; c < Para.NumberOfColours; c++)
        Transpose<1>(Data.col+c*size, Data.colslice+c*size, idim, jdim, kdim);
  }  
}

//
// Save subgrid flux on GPU and transfer the results back to CPU
//
extern "C"
void cuPPMSaveSubgridFluxes(
  fluxes *SubgridFluxes[], int NumberOfSubgrids, Elong_int GridGlobalStart[],
  cuPPMData &Data, cuPPMParameter &Para, int dir)
{
  int dim, idim, jdim;
  if (dir == 0) dim = 0, idim = 1, jdim = 2;
  else if (dir == 1) dim = 1, idim = 0, jdim = 2;
  else if (dir == 2) dim = 2, idim = 0, jdim = 1;
  const int size = 
    Para.GridDimension[0]*Para.GridDimension[1]*Para.GridDimension[2];
  int *FluxId = (int*)malloc(Para.MaxNumberOfBaryonFields*sizeof(int));
  float **Flux3D = (float**)malloc(Para.MaxNumberOfBaryonFields*sizeof(float*));
  int FluxCount = 5;
  FluxId[0] = Para.DensNum, FluxId[1] = Para.TENum;
  FluxId[2] = Para.Vel1Num, FluxId[3] = Para.Vel2Num, FluxId[4] = Para.Vel3Num;
  Flux3D[0] = Data.df, Flux3D[1] = Data.ef;
  if (dir == 0)
    Flux3D[2] = Data.uf, Flux3D[3] = Data.vf, Flux3D[4] = Data.wf;
  else if (dir == 1) 
    Flux3D[2] = Data.wf, Flux3D[3] = Data.uf, Flux3D[4] = Data.vf;
  else if (dir == 2) 
    Flux3D[2] = Data.vf, Flux3D[3] = Data.wf, Flux3D[4] = Data.uf;
  if (Para.DualEnergyFormalism) {
    FluxId[5] = Para.GENum;
    Flux3D[5] = Data.gef;
    FluxCount++;
  }
  for (int c = 0; c < Para.NumberOfColours; c++) {
    FluxId[6+c] = Para.colnum[c];
    Flux3D[6+c] = Data.colf + c*size;
    FluxCount++;
  }
  for (int n = 0; n < NumberOfSubgrids; n++) {
    const int fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] -
      GridGlobalStart[idim];
    const int fiend   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
      GridGlobalStart[idim];
    const int fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
      GridGlobalStart[jdim];
    const int fjend   = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
      GridGlobalStart[jdim];
    const int lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] - 
      GridGlobalStart[dim];
    const int rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
      GridGlobalStart[dim] + 1;
    int nfi = fiend - fistart + 1;
    int nfj = fjend - fjstart + 1;
    int nf = nfi*nfj;
    float *LeftFlux, *RightFlux;
    cumalloc((void**)&LeftFlux, nf*sizeof(float));
    cumalloc((void**)&RightFlux, nf*sizeof(float));

    for (int i = 0; i < FluxCount; i++)
      SaveSubgridFluxCUDA(SubgridFluxes[n]->LeftFluxes[FluxId[i]][dim],
                          SubgridFluxes[n]->RightFluxes[FluxId[i]][dim],
                          Flux3D[i], LeftFlux, RightFlux,
                          Para.GridDimension[0], 
                          Para.GridDimension[1], 
                          Para.GridDimension[2],
                          fistart, fiend, fjstart, fjend, lface, rface, dir);
    cufree(LeftFlux);
    cufree(RightFlux);
  }
  free(FluxId);
  free(Flux3D);
}




