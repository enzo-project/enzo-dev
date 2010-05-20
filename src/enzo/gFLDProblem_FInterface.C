/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class Fortran 
/  interfaces.
/
/  written by: Daniel Reynolds
/  date:       November, 2006
/  modified1:  July 20, 2007 by John Hayes; added Model to gfldproblem_matrixentries
/              and gfldproblem_diffrhs prototypes.
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(gfldproblem_matrixentries_3d)(
   Eflt64 *matentries, float *Er, float *Er0, float *Temp, float *kappaE, 
   float *adjvec, int *LimType, float *dt, FLOAT *a, float *theta, 
   float *aUnits, float *LenUnits, float *ErUnits, float *dx, 
   float *dy, float *dz, int *x0s, int *x0e, int *x1s, int *x1e, int *x2s, 
   int *x2e, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, int *xrface, int *yrface, int *zrface, 
   int *Model, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_matrixentries_2d)(
   Eflt64 *matentries, float *Er, float *Er0, float *Temp, float *kappaE, 
   float *adjvec, int *LimType, float *dt, FLOAT *a, float *theta, 
   float *aUnits, float *LenUnits, float *ErUnits, float *dx, 
   float *dy, int *x0s, int *x0e, int *x1s, int *x1e, int *Nx, int *Ny, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *xrface, int *yrface, 
   int *Model, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_matrixentries_1d)(
   Eflt64 *matentries, float *Er, float *Er0, float *Temp, float *kappaE, 
   float *adjvec, int *LimType, float *dt, FLOAT *a, float *theta, 
   float *aUnits, float *LenUnits, float *ErUnits, float *dx, 
   int *x0s, int *x0e, int *Nx, int *NGxl, int *NGxr, int *xrface, 
   int *Model, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_setnewtonbcs_3d)(
   Eflt64 *matentries, float *rhsentries, FLOAT *a, float *aUnits, 
   float *LenUnits, float *ErUnits, float *dx, float *dy, float *dz, 
   int *x0s, int *x0e, int *x1s, int *x1e, int *x2s, int *x2e, int *Nx, 
   int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, int *BCxL, int *BCxR, int *BCyL, int *BCyR, 
   int *BCzL, int *BCzR, int *xlface, int *xrface, int *ylface, 
   int *yrface, int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_setnewtonbcs_2d)(
   Eflt64 *matentries, float *rhsentries, FLOAT *a, float *aUnits, 
   float *LenUnits, float *ErUnits, float *dx, float *dy, int *x0s, 
   int *x0e, int *x1s, int *x1e, int *Nx, int *Ny, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *BCxL, int *BCxR, int *BCyL, 
   int *BCyR, int *xlface, int *xrface, int *ylface, int *yrface, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_setnewtonbcs_1d)(
   Eflt64 *matentries, float *rhsentries, FLOAT *a, float *aUnits, 
   float *LenUnits, float *ErUnits, float *dx, int *x0s, int *x0e, 
   int *Nx, int *NGxl, int *NGxr, int *BCxL, int *BCxR, int *xlface, 
   int *xrface, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_diffrhs_3d)(
   float *drhs, float *Er, float *Er0, float *Temp, float *kappaE, int *LimType, 
   FLOAT *a, float *aUnits, float *LenUnits, float *ErUnits, float *dx, 
   float *dy, float *dz, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *xlface, int *xrface, 
   int *ylface, int *yrface, int *zlface, int *zrface, int *Model, int *ier);
  
extern "C" void FORTRAN_NAME(gfldproblem_diffrhs_2d)(
   float *drhs, float *Er, float *Er0, float *Temp, float *kappaE, int *LimType, 
   FLOAT *a, float *aUnits, float *LenUnits, float *ErUnits, float *dx, 
   float *dy, int *Nx, int *Ny, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *xlface, int *xrface, int *ylface, int *yrface, int *Model, int *ier);
  
extern "C" void FORTRAN_NAME(gfldproblem_diffrhs_1d)(
   float *drhs, float *Er, float *Er0, float *Temp, float *kappaE, int *LimType, 
   FLOAT *a, float *aUnits, float *LenUnits, float *ErUnits, float *dx, int *Nx, 
   int *NGxl, int *NGxr, int *xlface, int *xrface, int *Model, int *ier);
  
extern "C" void FORTRAN_NAME(gfldproblem_localrhs)(
   float *Errhs, float *ecrhs, float *HIrhs, float *HeIrhs, float *HeIIrhs, 
   float *Ersrc, float *ecsrc, float *HIsrc, float *HeIsrc, float *HeIIsrc, 
   float *time, float *vx, float *vy, float *vz, float *rho, float *ec, 
   float *Er, float *nHI, float *nHeI, float *nHeII, float *Temp, float *eh, 
   float *kappaP, float *kappaE, FLOAT *a, FLOAT *adot, float *gamma, 
   float *hfrac, int *model, int *AnalyticChem, int *ESpectrum, float *CompA, 
   float *Comp_xray, float *Comp_temp, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, float *IsEsHeIInu, 
   int *NTempBins, float *TempStart, float *TempEnd, float *k1Tb, float *k2Tb, 
   float *k3Tb, float *k4Tb, float *k5Tb, float *k6Tb, float *ceHITb, 
   float *ceHeITb, float *ceHeIITb, float *ciHITb, float *ciHeITb, 
   float *ciHeISTb, float *ciHeIITb, float *reHIITb, float *reHeII1Tb, 
   float *reHeII2Tb, float *reHeIIITb, float *bremTb, float *piHI, float *piHeI, 
   float *piHeII, float *aUnits, float *DenUnits, float *VelUnits, 
   float *LenUnits, float *ErUnits, float *ecUnits, float *NiUnits, 
   float *ecScale, int *Nchem, float *dx, float *dy, float *dz, int *Nx, 
   int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, 
   int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_localjac)(
   float *Erjac_Er, float *Erjac_ec, float *Erjac_HI, float *Erjac_HeI, 
   float *Erjac_HeII, float *ecjac_Er, float *ecjac_ec, float *ecjac_HI, 
   float *ecjac_HeI, float *ecjac_HeII, float *HIjac_Er, float *HIjac_ec, 
   float *HIjac_HI, float *HIjac_HeI, float *HIjac_HeII, float *HeIjac_Er, 
   float *HeIjac_ec, float *HeIjac_HI, float *HeIjac_HeI, float *HeIjac_HeII, 
   float *HeIIjac_Er, float *HeIIjac_ec, float *HeIIjac_HI, float *HeIIjac_HeI, 
   float *HeIIjac_HeII, float *time, float *Er, float *ec, float *nHI, 
   float *nHeI, float *nHeII, float *eh, float *rho, float *vx, float *vy, 
   float *vz, int *Nchem, float *hfrac, int *model, int *ESpectrum, int *probtype, 
   int *dualenergy, FLOAT *a, FLOAT *adot, float *CompA, float *Comp_xray, 
   float *Comp_temp, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, float *IsEsHeIInu, 
   float *PlC0, float *PlC1, float *PlC2, float *PlC3, float *PlC4, 
   float *RoC0, float *RoC1, float *RoC2, float *RoC3, float *RoC4, 
   float *gamma, int *NTempBins, float *TempStart, float *TempEnd, 
   float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, 
   float *k6Tb, float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, 
   float *ciHeITb, float *ciHeISTb, float *ciHeIITb, float *reHIITb, 
   float *reHeII1Tb, float *reHeII2Tb, float *reHeIIITb, float *bremTb, 
   float *piHI, float *piHeI, float *piHeII, float *aUnits, float *DenUnits, 
   float *VelUnits, float *LenUnits, float *ErUnits, float *ecUnits, 
   float *NiUnits, 
   float *ecScale, float *dx, float *dy, float *dz, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(blocksolve)(
   float *Amat, float *xvec, float *bvec, int *N, int *M, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_opacity)(
   float *kappaP, float *kappaE, float *time, float *rho, 
   float *n_HI, float *n_HeI, float *n_HeII, float *Temperature, FLOAT *a, 
   int *Model, float *IsE, float *IsEsHI, float *IsEsHInu, float *IsEsHeI, 
   float *IsEsHeInu, float *IsEsHeII, float *IsEsHeIInu, float *PlC0, 
   float *PlC1, float *PlC2, float *PlC3, float *PlC4, float *EmC0, 
   float *EmC1, float *EmC2, float *EmC3, float *EmC4, float *aUnits, 
   float *DenUnits, float *LenUnits, float *TimeUnits, float *NiUnits, 
   int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, 
   float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_radiationsource)(
   float *Ersrc, float *time, float *Er, float *ec, float *n_HI, 
   float *n_HeI, float *n_HeII, float *Temperature, float *rho, float *eh, 
   float *vx, float *vy, float *vz, FLOAT *a, int *Model, 
   int *ProblemType, int *Nchem, float *HFrac, int *ESpectrum, 
   float *NGammaDot, float *EtaRadius, float *EtaCenter0, 
   float *EtaCenter1, float *EtaCenter2, float *aUnits, float *DenUnits, 
   float *VelUnits, float *LenUnits, float *TimeUnits, float *ErUnits, 
   float *ecUnits, float *NiUnits, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_gasenergysource)(
   float *ecsrc, float *time, float *Er, float *ec, float *n_HI, 
   float *n_HeI, float *n_HeII, float *Temperature, float *rho, float *eh, 
   float *vx, float *vy, float *vz, FLOAT *a, int *Model, int *ProblemType, 
   int *Nchem, float *HFrac, float *aUnits, float *DenUnits, float *VelUnits, 
   float *LenUnits, float *TimeUnits, float *ErUnits, float *ecUnits, 
   float *NiUnits, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_chemistrysource)(
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *time, float *Er, 
   float *ec, float *n_HI, float *n_HeI, float *n_HeII, float *Temperature, 
   float *rho, float *eh, float *vx, float *vy, float *vz, FLOAT *a, 
   int *Model, int *Nchem, float *HFrac, float *aUnits, float *DenUnits, 
   float *VelUnits, float *LenUnits, float *TimeUnits, float *ErUnits, 
   float *ecUnits, float *NiUnits, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_analyticinitguess)(
   float *Er, float *ec, float *HI, float *HeI, float *HeII, float *dt, float *vx, 
   float *vy, float *vz, float *rho, float *eh, float *src_Er, float *src_ec, 
   float *src_HI, float *src_HeI, float *src_HeII, float *gamma, float *HFrac, 
   int *Model, int *ESpectrum, int *DualEnergy, FLOAT *a, FLOAT *adot, float *CompA, 
   float *Comp_xray, float *Comp_temp, float *IsE, float *IsEsHI, float *IsEsHInu, 
   float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, float *IsEsHeIInu, 
   int *NTempBins, float *TempStart, float *TempEnd, float *k1Tb, float *k2Tb, 
   float *k3Tb, float *k4Tb, float *k5Tb, float *k6Tb, float *ceHITb, 
   float *ceHeITb, float *ceHeIITb, float *ciHITb, float *ciHeITb, float *ciHeISTb, 
   float *ciHeIITb, float *reHIITb, float *reHeII1Tb, float *reHeII2Tb, 
   float *reHeIIITb, float *bremTb, float *aUnits, float *DenUnits, 
   float *VelUnits, float *LenUnits, float *TimeUnits, float *ErUnits, 
   float *ecUnits, float *NiUnits, float *ecScale, int *Nchem, float *dx, 
   float *dy, float *dz, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);

extern "C" void FORTRAN_NAME(gfldproblem_analyticresid)(
   float *ecres, float *HIres, float *HeIres, float *HeIIres, float *Er, 
   float *ec, float *HI, float *HeI, float *HeII, float *Er0, float *ec0, 
   float *HI0, float *HeI0, float *HeII0, float *dt, float *vx, 
   float *vy, float *vz, float *rho, float *eh, float *src_Er, float *src_ec, 
   float *src_HI, float *src_HeI, float *src_HeII, float *gamma, float *HFrac, 
   int *Model, int *DualEnergy, FLOAT *a, FLOAT *adot, float *CompA, 
   float *Comp_xray, float *Comp_temp, float *IsE, float *IsEsHI, 
   float *IsEsHInu, float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, 
   float *IsEsHeIInu, int *NTempBins, float *TempStart, float *TempEnd, 
   float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, float *k6Tb, 
   float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, float *ciHeITb, 
   float *ciHeISTb, float *ciHeIITb, float *reHIITb, float *reHeII1Tb, 
   float *reHeII2Tb, float *reHeIIITb, float *bremTb, float *aUnits, float *DenUnits, 
   float *VelUnits, float *LenUnits, float *TimeUnits, float *ErUnits, 
   float *ecUnits, float *NiUnits, float *ecScale, int *Nchem, int *Nx, int *Ny, 
   int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *ier);



/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the gFLDProblem class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int gFLDProblem::MatrixEntries(Eflt64 *matentries, float *Er, 
			       float *Er0, float *Temperature, 
			       float *sigA, float *adjvec) 
{
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int x0s=1, x0e=LocDims[0], x1s=1, x1e=LocDims[1], x2s=1, x2e=LocDims[2];
  x0s -= (OnBdry[0][0] && (BdryType[0][0]==1)) ? 1 : 0;
  x0e += (OnBdry[0][1] && (BdryType[0][1]==1)) ? 1 : 0;
  x1s -= (OnBdry[1][0] && (BdryType[1][0]==1)) ? 1 : 0;
  x1e += (OnBdry[1][1] && (BdryType[1][1]==1)) ? 1 : 0;
  x2s -= (OnBdry[2][0] && (BdryType[2][0]==1)) ? 1 : 0;
  x2e += (OnBdry[2][1] && (BdryType[2][1]==1)) ? 1 : 0;
  int ier;
  if (rank == 3) {
    FORTRAN_NAME(gfldproblem_matrixentries_3d)
      (matentries, Er, Er0, Temperature, sigA, adjvec, &LimType, &dt, &a, 
       &theta, &aUnits, &LenUnits, &ErUnits, &dx[0], &dx[1], &dx[2], 
       &x0s, &x0e, 
       &x1s, &x1e, &x2s, &x2e, &LocDims[0], &LocDims[1], &LocDims[2], 
       &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
       &GhDims[2][0], &GhDims[2][1], &xrface, &yrface, &zrface, &Model, &ier);
    return(ier);
  }
  else if (rank == 2) {
    FORTRAN_NAME(gfldproblem_matrixentries_2d)
      (matentries, Er, Er0, Temperature, sigA, adjvec, &LimType, &dt, &a, 
       &theta, &aUnits, &LenUnits, &ErUnits, &dx[0], &dx[1], &x0s, 
       &x0e, &x1s, 
       &x1e, &LocDims[0], &LocDims[1], &GhDims[0][0], &GhDims[0][1], 
       &GhDims[1][0], &GhDims[1][1], &xrface, &yrface, &Model, &ier);
    return(ier);
  }
  else {
    FORTRAN_NAME(gfldproblem_matrixentries_1d)
      (matentries, Er, Er0, Temperature, sigA, adjvec, &LimType, &dt, &a, 
       &theta, &aUnits, &LenUnits, &ErUnits, &dx[0], &x0s, &x0e, 
       &LocDims[0], 
       &GhDims[0][0], &GhDims[0][1], &xrface, &Model, &ier);
    return(ier);
  }
}

/********/
int gFLDProblem::SetNewtonBCs(Eflt64 *matentries, float *rhsentries) 
{
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int x0s=1, x0e=LocDims[0], x1s=1, x1e=LocDims[1], x2s=1, x2e=LocDims[2];
  x0s -= (OnBdry[0][0] && (BdryType[0][0]==1)) ? 1 : 0;
  x0e += (OnBdry[0][1] && (BdryType[0][1]==1)) ? 1 : 0;
  x1s -= (OnBdry[1][0] && (BdryType[1][0]==1)) ? 1 : 0;
  x1e += (OnBdry[1][1] && (BdryType[1][1]==1)) ? 1 : 0;
  x2s -= (OnBdry[2][0] && (BdryType[2][0]==1)) ? 1 : 0;
  x2e += (OnBdry[2][1] && (BdryType[2][1]==1)) ? 1 : 0;
  int ier;
  if (rank == 3) {
    FORTRAN_NAME(gfldproblem_setnewtonbcs_3d)
      (matentries, rhsentries, &a, &aUnits, &LenUnits, &ErUnits, 
       &dx[0], &dx[1], &dx[2], &x0s, &x0e, &x1s, &x1e, &x2s, 
       &x2e, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
       &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
       &GhDims[2][1], &BdryType[0][0], &BdryType[0][1], &BdryType[1][0], 
       &BdryType[1][1], &BdryType[2][0], &BdryType[2][1], &xlface, 
       &xrface, &ylface, &yrface, &zlface, &zrface, &ier);
    return(ier);
  }
  else if (rank == 2) {
    FORTRAN_NAME(gfldproblem_setnewtonbcs_2d)
      (matentries, rhsentries, &a, &aUnits, &LenUnits, &ErUnits, 
       &dx[0], &dx[1], &x0s, &x0e, &x1s, &x1e, &LocDims[0], 
       &LocDims[1], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
       &GhDims[1][1], &BdryType[0][0], &BdryType[0][1], 
       &BdryType[1][0], &BdryType[1][1], &xlface, &xrface, &ylface, 
       &yrface, &ier);
    return(ier);
  }
  else {
    FORTRAN_NAME(gfldproblem_setnewtonbcs_1d)
      (matentries, rhsentries, &a, &aUnits, &LenUnits, &ErUnits, 
       &dx[0], &x0s, &x0e, &LocDims[0], &GhDims[0][0], &GhDims[0][1], 
       &BdryType[0][0], &BdryType[0][1], &xlface, &xrface, &ier);
    return(ier);
  }
}

/********/
int gFLDProblem::DiffRHS(float *drhs, float *Er, float *Er0, 
			 float *Temperature, float *sigA)
{
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int ier;
  if (rank == 3) {
    FORTRAN_NAME(gfldproblem_diffrhs_3d)
      (drhs, Er, Er0, Temperature, sigA, &LimType, &a, &aUnits, &LenUnits, 
       &ErUnits, &dx[0], &dx[1], &dx[2], &LocDims[0], &LocDims[1], 
       &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
       &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &xlface, &xrface, 
       &ylface, &yrface, &zlface, &zrface, &Model, &ier);
    return(ier);
  }
  else if (rank == 2) {
    FORTRAN_NAME(gfldproblem_diffrhs_2d)
      (drhs, Er, Er0, Temperature, sigA, &LimType, &a, &aUnits, &LenUnits, 
       &ErUnits, &dx[0], &dx[1], &LocDims[0], &LocDims[1], &GhDims[0][0], 
       &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &xlface, &xrface, 
       &ylface, &yrface, &Model, &ier);
    return(ier);
  }
  else {
    FORTRAN_NAME(gfldproblem_diffrhs_1d)
      (drhs, Er, Er0, Temperature, sigA, &LimType, &a, &aUnits, &LenUnits, 
       &ErUnits, &dx[0], &LocDims[0], &GhDims[0][0], &GhDims[0][1], 
       &xlface, &xrface, &Model, &ier);
    return(ier);
  }
}

/********/
int gFLDProblem::LocalRHS(float *rhs_Er, float *rhs_ec, float *rhs_HI, 
			  float *rhs_HeI, float *rhs_HeII, float *src_Er, 
			  float *src_ec, float *src_HI, float *src_HeI, 
			  float *src_HeII, float *time, float *ec, float *Er, 
			  float *Temperature, float *kappaP, float *kappaE, 
			  float *nHI, float *nHeI, float *nHeII)
{
  int ier;
  FORTRAN_NAME(gfldproblem_localrhs)
    (rhs_Er, rhs_ec, rhs_HI, rhs_HeI, rhs_HeII, src_Er, src_ec, src_HI, src_HeI, 
     src_HeII, time, vx, vy, vz, rho, ec, Er, nHI, nHeI, nHeII, Temperature, eh, 
     kappaP, kappaE, &a, &adot, &Gamma, &HFrac, &Model, &AnalyticChem, &ESpectrum, 
     &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, &intSigE, 
     &intSigESigHI, &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, 
     &intSigESigHeII, &intSigESigHeIInu, &CoolData.NumberOfTemperatureBins, 
     &CoolData.TemperatureStart, &CoolData.TemperatureEnd, RateData.k1, 
     RateData.k2, RateData.k3, RateData.k4, RateData.k5, RateData.k6, 
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI, 
     CoolData.ciHeI, CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
     CoolData.reHeII1, CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, 
     &CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII, &aUnits, &DenUnits, 
     &VelUnits, &LenUnits, &ErUnits, &ecUnits, &NiUnits, &ecScale, &Nchem, 
     &dx[0], &dx[1], 
     &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::LocalJac(float *Erjac_Er, float *Erjac_ec, float *Erjac_HI, 
			  float *Erjac_HeI, float *Erjac_HeII, float *ecjac_Er, 
			  float *ecjac_ec, float *ecjac_HI, float *ecjac_HeI, 
			  float *ecjac_HeII, float *HIjac_Er, float *HIjac_ec, 
			  float *HIjac_HI, float *HIjac_HeI, float *HIjac_HeII, 
			  float *HeIjac_Er, float *HeIjac_ec, float *HeIjac_HI, 
			  float *HeIjac_HeI, float *HeIjac_HeII, float *HeIIjac_Er, 
			  float *HeIIjac_ec, float *HeIIjac_HI, float *HeIIjac_HeI, 
			  float *HeIIjac_HeII, float *time, float *Er, float *ec, 
			  float *nHI, float *nHeI, float *nHeII)
{
  int ier;
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  FORTRAN_NAME(gfldproblem_localjac)
    (Erjac_Er, Erjac_ec, Erjac_HI, Erjac_HeI, Erjac_HeII, ecjac_Er, ecjac_ec, 
     ecjac_HI, ecjac_HeI, ecjac_HeII, HIjac_Er, HIjac_ec, HIjac_HI, HIjac_HeI, 
     HIjac_HeII, HeIjac_Er, HeIjac_ec, HeIjac_HI, HeIjac_HeI, HeIjac_HeII, 
     HeIIjac_Er, HeIIjac_ec, HeIIjac_HI, HeIIjac_HeI, HeIIjac_HeII, time, Er, 
     ec, nHI, nHeI, nHeII, eh, rho, vx, vy, vz, &Nchem, &HFrac, &Model, 
     &ESpectrum, &ProblemType, &dualenergy, &a, &adot, &CoolData.comp, 
     &CoolData.comp_xray, &CoolData.temp_xray, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, &intSigESigHeII, 
     &intSigESigHeIInu, &PlanckOpacityC0, &PlanckOpacityC1, &PlanckOpacityC2, 
     &PlanckOpacityC3, &PlanckOpacityC4, &EnergyOpacityC0, &EnergyOpacityC1, 
     &EnergyOpacityC2, &EnergyOpacityC3, &EnergyOpacityC4, &Gamma, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &CoolData.piHI, &CoolData.piHeI, 
     &CoolData.piHeII, &aUnits, &DenUnits, &VelUnits, &LenUnits, &ErUnits, 
     &ecUnits, 
     &NiUnits, &ecScale, &dx[0], &dx[1], &dx[2], &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::BlockSolve(float *Amat, float *xvec, 
			    float *bvec, int *N, int *M)
{
  int ier;
  FORTRAN_NAME(blocksolve)(Amat, xvec, bvec, N, M, &ier);
  return(ier);
}


/********/
int gFLDProblem::Opacity(float *kappaP, float *kappaE, float *time, 
			 float *n_HI, float *n_HeI, float *n_HeII, 
			 float *Temperature)
{
  int ier;
  FORTRAN_NAME(gfldproblem_opacity)
    (kappaP, kappaE, time, rho, n_HI, n_HeI, n_HeII, 
     Temperature, &a, &Model, &intSigE, &intSigESigHI, &intSigESigHInu, 
     &intSigESigHeI, &intSigESigHeInu, &intSigESigHeII, &intSigESigHeIInu, 
     &PlanckOpacityC0, &PlanckOpacityC1, &PlanckOpacityC2, 
     &PlanckOpacityC3, &PlanckOpacityC4, &EnergyOpacityC0, &EnergyOpacityC1, 
     &EnergyOpacityC2, &EnergyOpacityC3, &EnergyOpacityC4, 
     &aUnits, &DenUnits, &LenUnits, &TimeUnits, &NiUnits, &Nchem, &LocDims[0], 
     &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::RadiationSource(float *Ersrc, float *time, float *Er, 
				 float *ec, float *n_HI, float *n_HeI, 
				 float *n_HeII, float *Temperature)
{
  int ier;
  FORTRAN_NAME(gfldproblem_radiationsource)
    (Ersrc, time, Er, ec, n_HI, n_HeI, n_HeII, Temperature, rho, eh, vx, vy, vz, 
     &a, &Model, &ProblemType, &Nchem, &HFrac, &ESpectrum, &IonizationParms[0], 
     &IonizationParms[1], &IonizationParms[2], &IonizationParms[3], &IonizationParms[4], 
     &aUnits, &DenUnits, &VelUnits, &LenUnits, &TimeUnits, &ErUnits, &ecUnits, &NiUnits, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], 
     &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::GasEnergySource(float *ecsrc, float *time, float *Er, 
				 float *ec, float *n_HI, float *n_HeI, 
				 float *n_HeII, float *Temperature)
{
  int ier;
  FORTRAN_NAME(gfldproblem_gasenergysource)
    (ecsrc, time, Er, ec, n_HI, n_HeI, n_HeII, Temperature, rho, eh, vx, vy, vz, &a, 
     &Model, &ProblemType, &Nchem, &HFrac, &aUnits, &DenUnits, &VelUnits, &LenUnits, 
     &TimeUnits, &ErUnits, &ecUnits, &NiUnits, &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::ChemistrySource(float *HIsrc, float *HeIsrc, float *HeIIsrc, 
				 float *time, float *Er, float *ec, float *n_HI, 
				 float *n_HeI, float *n_HeII, float *Temperature)
{
  int ier;
  FORTRAN_NAME(gfldproblem_chemistrysource)
    (HIsrc, HeIsrc, HeIIsrc, time, Er, ec, n_HI, n_HeI, n_HeII, Temperature, 
     rho, eh, vx, vy, vz, &a, &Model, &Nchem, &HFrac, 
     &aUnits, &DenUnits, &VelUnits, &LenUnits, &TimeUnits, &ErUnits, &ecUnits, &NiUnits, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], 
     &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::AnalyticInitGuess(EnzoVector *u, float deltat)
{
  int ier;
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  float *Er   = u->GetData(0);
  float *ec   = u->GetData(1);
  float *HI   = u->GetData(2);
  float *HeI  = u->GetData(3);
  float *HeII = u->GetData(4);
  float *Ersrc   = extsrc->GetData(0);
  float *ecsrc   = extsrc->GetData(1);
  float *HIsrc   = extsrc->GetData(2);
  float *HeIsrc  = extsrc->GetData(3);
  float *HeIIsrc = extsrc->GetData(4);
  FORTRAN_NAME(gfldproblem_analyticinitguess)
    (Er, ec, HI, HeI, HeII, &deltat, vx, vy, vz, rho, eh, Ersrc, ecsrc, 
     HIsrc, HeIsrc, HeIIsrc, &Gamma, &HFrac, &Model, &ESpectrum, &dualenergy, 
     &a, &adot, &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &intSigE, &intSigESigHI, &intSigESigHInu, &intSigESigHeI, 
     &intSigESigHeInu, &intSigESigHeII, &intSigESigHeIInu, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &aUnits, &DenUnits, &VelUnits, 
     &LenUnits, &TimeUnits, &ErUnits, &ecUnits, &NiUnits, &ecScale, &Nchem, 
     &dx[0], &dx[1], &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDProblem::AnalyticResid(EnzoVector *u0, EnzoVector *u, 
			       EnzoVector *fu, float deltat)
{
  int ier;
  int dualenergy = (DualEnergyFormalism) ? 1 : 0;
  float *Er0     = u0->GetData(0);
  float *ec0     = u0->GetData(1);
  float *HI0     = u0->GetData(2);
  float *HeI0    = u0->GetData(3);
  float *HeII0   = u0->GetData(4);
  float *Er      = u->GetData(0);
  float *ec      = u->GetData(1);
  float *HI      = u->GetData(2);
  float *HeI     = u->GetData(3);
  float *HeII    = u->GetData(4);
  float *ecres   = fu->GetData(1);
  float *HIres   = fu->GetData(2);
  float *HeIres  = fu->GetData(3);
  float *HeIIres = fu->GetData(4);
  float *Ersrc   = extsrc->GetData(0);
  float *ecsrc   = extsrc->GetData(1);
  float *HIsrc   = extsrc->GetData(2);
  float *HeIsrc  = extsrc->GetData(3);
  float *HeIIsrc = extsrc->GetData(4);
  FORTRAN_NAME(gfldproblem_analyticresid)
    (ecres, HIres, HeIres, HeIIres, Er, ec, HI, HeI, HeII, Er0, ec0, HI0, 
     HeI0, HeII0, &deltat, vx, vy, vz, rho, eh, Ersrc, ecsrc, HIsrc, HeIsrc, 
     HeIIsrc, &Gamma, &HFrac, &Model, &dualenergy, &a, &adot, &CoolData.comp, 
     &CoolData.comp_xray, &CoolData.temp_xray, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, &intSigESigHeII, 
     &intSigESigHeIInu, &CoolData.NumberOfTemperatureBins, 
     &CoolData.TemperatureStart, &CoolData.TemperatureEnd, RateData.k1, 
     RateData.k2, RateData.k3, RateData.k4, RateData.k5, RateData.k6, 
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI, 
     CoolData.ciHeI, CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
     CoolData.reHeII1, CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, 
     &aUnits, &DenUnits, &VelUnits, &LenUnits, &TimeUnits, &ErUnits, 
     &ecUnits, &NiUnits, &ecScale, &Nchem, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], 
     &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}



#endif
