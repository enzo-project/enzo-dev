/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  Fortran interfaces
/
/  written by: Daniel Reynolds
/  date:       July 2009
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(gfldsplit_setupsystem)(
   Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, float *E0, float *E, float *Temp, 
   float *Temp0, float *kappaE, float *eta, float *dt, FLOAT *a, 
   FLOAT *a0, FLOAT *adot, FLOAT *adot0, int *ESpectrum, float *theta, 
   float *aUnits, float *LenUnits, float *LenUnits0, float *ErUnits, 
   float *ErUnits0, float *NiUnits, float *NiUnits0, int *rank, float *dx, 
   float *dy, float *dz, int *BCTypeXl, int *BCTypeXr, int *BCTypeYl, 
   int *BCTypeYr, int *BCTypeZl, int *BCTypeZr, int *x0s, int *x0e, int *x1s, 
   int *x1e, int *x2s, int *x2e, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, 
   int *NGyl, int *NGyr, int *NGzl, int *NGzr, int *xlface, int *xrface, 
   int *ylface, int *yrface, int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(gfldsplit_opacity)(
   float *kappaE, float *time, float *rho, float *n_HI, float *n_HeI, 
   float *n_HeII, FLOAT *a, int *Model, float *IsE, float *IsEsHI, 
   float *IsEsHInu, float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, 
   float *IsEsHeIInu, float *EmC0, float *EmC1, float *EmC2, float *aUnits, 
   float *DenUnits, float *LenUnits, float *TimeUnits, float *NiUnits, 
   int *Nchem, int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, 
   int *NGyr, int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, 
   float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldsplit_radiationsource)(
   float *Ersrc, float *time, FLOAT *a, int *ProblemType, int *ESpectrum, 
   float *NGammaDot, float *EtaRadius, float *EtaCenter, float *aUnits, 
   float *LenUnits, float *TimeUnits, float *ErUnits, int *Nx, int *Ny, int *Nz, 
   int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldsplit_gasenergysource)(
   float *ecsrc, float *time, FLOAT *a, int *ProblemType, float *aUnits, 
   float *VelUnits, float *LenUnits, float *TimeUnits, float *ecUnits, 
   int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, float *x0L, float *x0R, float *x1L, float *x1R, 
   float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldsplit_chemistrysource)(
   float *HIsrc, float *HeIsrc, float *HeIIsrc, float *time, FLOAT *a, 
   int *ProblemType, int *Nchem, float *HFrac, float *aUnits, float *LenUnits, 
   float *TimeUnits, float *NiUnits, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(gfldsplit_analyticinitguess)(
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

extern "C" void FORTRAN_NAME(gfldsplit_analyticchemistry)(
   float *Er, float *ec, float *HI, float *HeI, float *HeII, float *Er0, 
   float *ec0, float *HI0, float *HeI0, float *HeII0, float *dt, float *vx, 
   float *vy, float *vz, float *rho, float *eh, float *src_ec, float *src_HI, 
   float *src_HeI, float *src_HeII, float *kappa, float *gamma, float *HFrac, 
   int *Model, int *ProbType, int *DualEnergy, FLOAT *a, FLOAT *adot, 
   float *CompA, float *Comp_xray, float *Comp_temp, float *IsE, float *IsEsHI, 
   float *IsEsHInu, float *IsEsHeI, float *IsEsHeInu, float *IsEsHeII, 
   float *IsEsHeIInu, int *NTempBins, float *TempStart, float *TempEnd, 
   float *k1Tb, float *k2Tb, float *k3Tb, float *k4Tb, float *k5Tb, float *k6Tb, 
   float *ceHITb, float *ceHeITb, float *ceHeIITb, float *ciHITb, float *ciHeITb, 
   float *ciHeISTb, float *ciHeIITb, float *reHIITb, float *reHeII1Tb, 
   float *reHeII2Tb, float *reHeIIITb, float *bremTb, float *aUnits, 
   float *DenUnits, float *VelUnits, float *LenUnits, float *TimeUnits, 
   float *ErUnits, float *ecUnits, float *NiUnits, float *ecScale, int *Nchem, 
   int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, int *ier);



/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the gFLDSplit class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int gFLDSplit::SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, float *E0, 
			   float *E, float *kappaE, float *Temp, float *Temp0, 
			   float *eta) 
{
  int xlface = (OnBdry[0][0]) ? 1 : 0;
  int xrface = (OnBdry[0][1]) ? 1 : 0;
  int ylface = (OnBdry[1][0]) ? 1 : 0;
  int yrface = (OnBdry[1][1]) ? 1 : 0;
  int zlface = (OnBdry[2][0]) ? 1 : 0;
  int zrface = (OnBdry[2][1]) ? 1 : 0;
  int x0s=1, x0e=LocDims[0], x1s=1, x1e=LocDims[1], x2s=1, x2e=LocDims[2];
  int ier;
  FORTRAN_NAME(gfldsplit_setupsystem)
    (mat, rhs, rhsnorm, E0, E, Temp, Temp0, kappaE, eta, &dt, &a, &a0, 
     &adot, &adot0, &ESpectrum, &theta, &aUnits, &LenUnits, &LenUnits0, 
     &ErUnits, &ErUnits0, &NiUnits, &NiUnits0, &rank, &dx[0], &dx[1], 
     &dx[2], &BdryType[0][0], &BdryType[0][1], &BdryType[1][0], &BdryType[1][1], 
     &BdryType[2][0], &BdryType[2][1], &x0s, &x0e, &x1s, &x1e, &x2s, &x2e, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &xlface, 
     &xrface, &ylface, &yrface, &zlface, &zrface, &ier);

  // combine the processor-local rhsnorm values together before returning 
  // (since I perform no MPI in F90 modules)
  float rhssum=0.0;
#ifdef USE_MPI
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(rhsnorm, &rhssum, 1, DataType, MPI_SUM, MPI_COMM_WORLD);
#else
  rhssum = *rhsnorm;
#endif
  *rhsnorm = sqrt(rhssum);

  return(ier);
}


/********/
int gFLDSplit::Opacity(float *kappaE, float *time, EnzoVector *u)
{
  int ier;
  float *n_HI   = NULL;
  float *n_HeI  = NULL;
  float *n_HeII = NULL;
  if (Nchem > 0)  n_HI = u->GetData(2);
  if (Nchem > 1)  {
    n_HeI  = u->GetData(3);
    n_HeII = u->GetData(4);
  }
  FLOAT aval = (a+a0)*0.5;
  float dUn  = (DenUnits + DenUnits0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  float nUn  = (NiUnits + NiUnits0)*0.5;
  FORTRAN_NAME(gfldsplit_opacity)
    (kappaE, time, rho, n_HI, n_HeI, n_HeII, &aval, &Model, &intSigE, 
     &intSigESigHI, &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, 
     &intSigESigHeII, &intSigESigHeIInu, &EnergyOpacityC0, &EnergyOpacityC1, 
     &EnergyOpacityC2, &aUnits, &dUn, &lUn, &TimeUnits, &nUn, &Nchem, 
     &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], 
     &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], 
     &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDSplit::RadiationSource(float *Ersrc, float *time)
{
  int ier;
  FLOAT aval = (a+a0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  float rUn  = (ErUnits + ErUnits0)*0.5;
  FORTRAN_NAME(gfldsplit_radiationsource)
    (Ersrc, time, &aval, &ProblemType, &ESpectrum, &NGammaDot, &EtaRadius, 
     EtaCenter, &aUnits, &lUn, &TimeUnits, &rUn, &LocDims[0], &LocDims[1], 
     &LocDims[2], &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], 
     &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDSplit::GasEnergySource(float *ecsrc, float *time)
{
  int ier;
  FLOAT aval = (a+a0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  FORTRAN_NAME(gfldsplit_gasenergysource)
    (ecsrc, time, &aval, &ProblemType, &aUnits, &VelUnits, &lUn, 
     &TimeUnits, &ecUnits, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], 
     &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], &EdgeVals[1][1], 
     &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDSplit::ChemistrySource(float *HIsrc, float *HeIsrc, 
			       float *HeIIsrc, float *time)
{
  int ier;
  FLOAT aval = (a+a0)*0.5;
  float lUn  = (LenUnits + LenUnits0)*0.5;
  float nUn  = (NiUnits + NiUnits0)*0.5;
  FORTRAN_NAME(gfldsplit_chemistrysource)
    (HIsrc, HeIsrc, HeIIsrc, time, &aval, &ProblemType, &Nchem, &HFrac, 
     &aUnits, &lUn, &TimeUnits, &nUn, &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], &EdgeVals[1][0], 
     &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}


/********/
int gFLDSplit::AnalyticInitGuess(EnzoVector *u, float deltat)
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
  FLOAT aval = a;
  FLOAT adval = adot;
  float dUn  = DenUnits;
  float lUn  = LenUnits;
  float rUn  = ErUnits;
  float nUn  = NiUnits;
  FORTRAN_NAME(gfldsplit_analyticinitguess)
    (Er, ec, HI, HeI, HeII, &deltat, vx, vy, vz, rho, eh, Ersrc, ecsrc, 
     HIsrc, HeIsrc, HeIIsrc, &Gamma, &HFrac, &Model, &ESpectrum, &dualenergy, 
     &aval, &adval, &CoolData.comp, &CoolData.comp_xray, &CoolData.temp_xray, 
     &intSigE, &intSigESigHI, &intSigESigHInu, &intSigESigHeI, 
     &intSigESigHeInu, &intSigESigHeII, &intSigESigHeIInu, 
     &CoolData.NumberOfTemperatureBins, &CoolData.TemperatureStart, 
     &CoolData.TemperatureEnd, RateData.k1, RateData.k2, RateData.k3, 
     RateData.k4, RateData.k5, RateData.k6, CoolData.ceHI, CoolData.ceHeI, 
     CoolData.ceHeII, CoolData.ciHI, CoolData.ciHeI, CoolData.ciHeIS, 
     CoolData.ciHeII, CoolData.reHII, CoolData.reHeII1, CoolData.reHeII2, 
     CoolData.reHeIII, CoolData.brem, &aUnits, &dUn, &VelUnits, &lUn, 
     &TimeUnits, &rUn, &ecUnits, &nUn, &ecScale, &Nchem, &dx[0], &dx[1], 
     &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], 
     &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], 
     &GhDims[2][1], &ier);
  return(ier);
}


/********/
int gFLDSplit::AnalyticChemistry(EnzoVector *u0, EnzoVector *u, 
				 EnzoVector *src, float deltat)
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
  float *ecsrc   = src->GetData(1);
  float *HIsrc   = src->GetData(2);
  float *HeIsrc  = src->GetData(3);
  float *HeIIsrc = src->GetData(4);
  FLOAT aval = a;
  FLOAT adval = adot;
  float dUn  = DenUnits;
  float lUn  = LenUnits;
  float rUn  = ErUnits;
  float nUn  = NiUnits;
  FORTRAN_NAME(gfldsplit_analyticchemistry)
    (Er, ec, HI, HeI, HeII, Er0, ec0, HI0, HeI0, HeII0, &deltat, vx, vy, 
     vz, rho, eh, ecsrc, HIsrc, HeIsrc, HeIIsrc, OpacityE, &Gamma, 
     &HFrac, &Model, &ProblemType, &dualenergy, &aval, &adval, &CoolData.comp, 
     &CoolData.comp_xray, &CoolData.temp_xray, &intSigE, &intSigESigHI, 
     &intSigESigHInu, &intSigESigHeI, &intSigESigHeInu, &intSigESigHeII, 
     &intSigESigHeIInu, &CoolData.NumberOfTemperatureBins, 
     &CoolData.TemperatureStart, &CoolData.TemperatureEnd, RateData.k1, 
     RateData.k2, RateData.k3, RateData.k4, RateData.k5, RateData.k6, 
     CoolData.ceHI, CoolData.ceHeI, CoolData.ceHeII, CoolData.ciHI, 
     CoolData.ciHeI, CoolData.ciHeIS, CoolData.ciHeII, CoolData.reHII, 
     CoolData.reHeII1, CoolData.reHeII2, CoolData.reHeIII, CoolData.brem, 
     &aUnits, &dUn, &VelUnits, &lUn, &TimeUnits, &rUn, &ecUnits, &nUn, &ecScale, 
     &Nchem, &LocDims[0], &LocDims[1], &LocDims[2], &GhDims[0][0], &GhDims[0][1], 
     &GhDims[1][0], &GhDims[1][1], &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}

#endif
