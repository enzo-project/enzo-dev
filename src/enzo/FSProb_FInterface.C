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
/  Free-streaming Radiation Implicit Problem Class
/  Fortran interfaces.
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified :  
/
/  PURPOSE: provides C++ interfaces to the relevant Fortran 
/           computational kernels
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"

/* Fortran function prototypes */
extern "C" void FORTRAN_NAME(fsprob_setupsystem)(
   Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, float *E, float *E0,
   int *kappa_h2on, float *kappa_arr, float *kappa_c, float *eta, float *dt, 
   FLOAT *a, FLOAT *a0, FLOAT *adot, FLOAT *adot0, float *theta, 
   float *LenUnits, float *LenUnits0, float *EUnits, float *EUnits0, 
   float *DUnits, float *Dunits0, int *rank, float *dx, float *dy, float *dz, 
   int *BCTypeXl, int *BCTypeXr, int *BCTypeYl, int *BCTypeYr, int *BCTypeZl, 
   int *BCTypeZr, int *x0s, int *x0e, int *x1s, int *x1e, int *x2s, int *x2e, 
   int *Nx, int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, 
   int *NGzl, int *NGzr, int *xlface, int *xrface, int *ylface, int *yrface, 
   int *zlface, int *zrface, int *ier);

extern "C" void FORTRAN_NAME(fsprob_radiationsource)(
   float *Efsrc, float *time, FLOAT *a, int *ProblemType, Eflt64 *NGammaDot, 
   float *EtaRadius, float *EtaCenter, float *aUnits, float *LenUnits, 
   float *TimeUnits, float *EUnits, int *Nx, int *Ny, int *Nz, int *NGxl, 
   int *NGxr, int *NGyl, int *NGyr, int *NGzl, int *NGzr, float *x0L, 
   float *x0R, float *x1L, float *x1R, float *x2L, float *x2R, int *ier);

extern "C" void FORTRAN_NAME(fsprob_initialguess)(
   float *Ef, float *Ef0, float *src_Ef, int *iguess, float *dt, 
   int *kappa_h2on, float *kappa, float *kappa_c,
   FLOAT *a, FLOAT *adot, float *aUnits, float *LenUnits, float *TimeUnits, 
   float *EUnits, float *DenUnits, float *dx, float *dy, float *dz, int *Nx, 
   int *Ny, int *Nz, int *NGxl, int *NGxr, int *NGyl, int *NGyr, int *NGzl, 
   int *NGzr, int *ier);




/* C++ interface wrappers */
/* These are each designed to extract as much relevant information 
   as possible from the FSProb class, so that the argument 
   lists are simplified to only the relevant input/output arguments. */

/********/
int FSProb::SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, 
			float *E, float *E0, float *eta, float *opacity)
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
  FORTRAN_NAME(fsprob_setupsystem)
    (mat, rhs, rhsnorm, E, E0, &kappa_h2on, opacity, &kappa0, eta, &dt, 
     &a, &a0, &adot, &adot0, &theta, &LenUnits, &LenUnits0, &EUnits, 
     &EUnits0, &DenUnits, &DenUnits0, &rank, &dx[0], &dx[1], &dx[2], 
     &BdryType[0][0], &BdryType[0][1], &BdryType[1][0], &BdryType[1][1], 
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
int FSProb::RadiationSource(float *Efsrc)
{
  int ier;
  float time = (told + tnew)/2.0;
  FLOAT aval = (a + a0)/2.0;
  float lUn  = (LenUnits + LenUnits0)/2.0;
  float tUn  = (TimeUnits + TimeUnits0)/2.0;
  float EUn  = (EUnits + EUnits0)/2.0;
  FORTRAN_NAME(fsprob_radiationsource)
    (Efsrc, &time, &aval, &ProblemType, &NGammaDot, &EtaRadius, EtaCenter, 
     &aUnits, &lUn, &tUn, &EUn, &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &EdgeVals[0][0], &EdgeVals[0][1], 
     &EdgeVals[1][0], &EdgeVals[1][1], &EdgeVals[2][0], &EdgeVals[2][1], &ier);
  return(ier);
}

/********/
int FSProb::InitialGuess(EnzoVector *Ef, EnzoVector *Ef0, EnzoVector *Efsrc)
{
  int ier;
  float *Efptr = Ef->GetData(0);
  float *Ef0ptr = Ef0->GetData(0);
  float *Efsrcptr = Efsrc->GetData(0);
  float *kappaptr = kappa->GetData(0);
  FORTRAN_NAME(fsprob_initialguess)
    (Efptr, Ef0ptr, Efsrcptr, &initial_guess, &dt, &kappa_h2on, kappaptr, 
     &kappa0, &a, &adot, &aUnits, &LenUnits, &TimeUnits, &EUnits, &DenUnits, 
     &dx[0], &dx[1], &dx[2], &LocDims[0], &LocDims[1], &LocDims[2], 
     &GhDims[0][0], &GhDims[0][1], &GhDims[1][0], &GhDims[1][1], 
     &GhDims[2][0], &GhDims[2][1], &ier);
  return(ier);
}

#endif
