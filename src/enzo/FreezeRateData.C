/***********************************************************************
/
/  INITIALIZE THE MULTI-SPECIES RATES
/
/  written by: Daniel R. Reynolds
/  date:       January 2011
/  modified1:  
/
/  PURPOSE:
/    For isothermal runs, this routine freezes all cooling/chemistry 
/    rates to the current temperature, through modifying the CoolData 
/    and RateData rate tables to hold a constant value.  This routine 
/    assumes that the user has a homogeneous temperature field, given
/    using the first active grid cell in the first local grid on this 
/    MPI process.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int InitializeRateData(FLOAT Time);



int FreezeRateData(FLOAT Time, HierarchyEntry &TopGrid)
{

  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.");

  // get current units
  double mUnit=1;
  float dUnit=1, lUnit=1, TUnit=1, tUnit=1, vUnit=1;
  if (GetUnits(&dUnit,&lUnit,&TUnit,&tUnit,&vUnit,&mUnit,Time) == FAIL) 
    ENZO_FAIL("Error in GetUnits.\n");
    
  // Find a TopLevel grid owned by this processor, to compute temperature
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) 
    ENZO_FAIL("FreezeRateData could not find a level 0 grid!");

  // compute temperature
  //    find the first active cell on this grid (idx)
  int rank = ThisGrid->GridData->GetGridRank();
  int ghZl = (rank > 2) ? NumberOfGhostZones : 0;
  int ghYl = (rank > 1) ? NumberOfGhostZones : 0;
  int ghXl = NumberOfGhostZones;
  int n3[] = {1, 1, 1};
  for (int dim=0; dim<rank; dim++)
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
  int idx = (ghZl*(n3[1] + 2*ghYl) + ghYl)*(n3[0] + 2*ghXl) + ghXl;

  //    compute the internal energy
  float eint;
  if (DualEnergyFormalism) {  // get internal energy
    float *eint_ = ThisGrid->GridData->AccessGasEnergy();
    if (eint_ == NULL)
      ENZO_FAIL("Could not access internal energy field!");
    eint = eint_[idx]*vUnit*vUnit;
  } else {   // must compute internal energy from total energy
    float *etot_ = ThisGrid->GridData->AccessTotalEnergy();
    if (etot_ == NULL)
      ENZO_FAIL("Could not access total energy field!");
    float *vx_ = ThisGrid->GridData->AccessVelocity1();
    float vx = (vx_ == NULL) ? 0.0 : vx_[idx];
    float *vy_ = ThisGrid->GridData->AccessVelocity2();
    float vy = (vy_ == NULL) ? 0.0 : vy_[idx];
    float *vz_ = ThisGrid->GridData->AccessVelocity3();
    float vz = (vz_ == NULL) ? 0.0 : vz_[idx];
    eint = vUnit*vUnit*(etot_[idx] - 0.5*(vx_[idx]*vx_[idx] 
			+ vy_[idx]*vy_[idx] + vz_[idx]*vz_[idx]));
  }

  //    compute the number density
  float *rho_ = ThisGrid->GridData->AccessDensity();
  if (rho_ == NULL)
    ENZO_FAIL("Could not access Density field!");
  float rho = rho_[idx];
  float *HI_ = ThisGrid->GridData->AccessHIDensity();
  if (HI_ == NULL)
    ENZO_FAIL("Could not access Density field!");
  float HI = HI_[idx];
  float *HeI_  = ThisGrid->GridData->AccessHeIDensity();
  float HeI = (HeI_ == NULL) ? 0.0 : HeI_[idx];
  float *HeII_ = ThisGrid->GridData->AccessHeIIDensity();
  float HeII = (HeII_ == NULL) ? 0.0 : HeII_[idx];
  float nHI = HI;
  float nHII = (rho*CoolData.HydrogenFractionByMass - nHI);
  float nHeI = HeI;
  float nHeII = HeII;
  float nHeIII = (rho*(1.0-CoolData.HydrogenFractionByMass) - nHeI - nHeII);
  float ne = nHII + nHeII/4.0 + nHeIII/2.0;
  float num_density = 0.25*(nHeI+nHeII+nHeIII) + nHI + nHII + ne;

  //    put the temperature together
  float mu = rho/num_density;
  float mp = 1.67262171e-24;    // proton mass [g]
  float kb = 1.3806504e-16;     // boltzmann constant [erg/K]
  float Temp = max((Gamma-1.0)*mu*mp*eint/kb, 1.0);

  // find temperature bin
  float lamT = 3.15614e5/Temp;
  float lTempS = log(CoolData.TemperatureStart);
  float lTempE = log(CoolData.TemperatureEnd);
  float dlTemp = (lTempE - lTempS)/(1.0*CoolData.NumberOfTemperatureBins - 1.0);
  float lTemp  = min(max(log(Temp), lTempS), lTempE);
  int Tidx = min(CoolData.NumberOfTemperatureBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1));
  int Tidxp = Tidx+1;
  float Tl = lTempS + (Tidx-1)*dlTemp;
  float Tr = lTempS +  Tidx*dlTemp;
  float Tfac = (lTemp - Tl)/(Tr - Tl);


  // update all rate tables to use fixed values for this temperature
  float rate = CoolData.ceHI[Tidx] + (CoolData.ceHI[Tidxp] - CoolData.ceHI[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ceHI[i] = rate;
  
  rate = CoolData.ceHeI[Tidx] + (CoolData.ceHeI[Tidxp] - CoolData.ceHeI[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ceHeI[i] = rate;
  
  rate = CoolData.ceHeII[Tidx] + (CoolData.ceHeII[Tidxp] - CoolData.ceHeII[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ceHeII[i] = rate;
  
  rate = CoolData.ciHI[Tidx] + (CoolData.ciHI[Tidxp] - CoolData.ciHI[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ciHI[i] = rate;
  
  rate = CoolData.ciHeI[Tidx] + (CoolData.ciHeI[Tidxp] - CoolData.ciHeI[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ciHeI[i] = rate;
  
  rate = CoolData.ciHeIS[Tidx] + (CoolData.ciHeIS[Tidxp] - CoolData.ciHeIS[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ciHeIS[i] = rate;
  
  rate = CoolData.ciHeII[Tidx] + (CoolData.ciHeII[Tidxp] - CoolData.ciHeII[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.ciHeII[i] = rate;
  
  rate = CoolData.reHII[Tidx] + (CoolData.reHII[Tidxp] - CoolData.reHII[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.reHII[i] = rate;
  
  rate = CoolData.reHeII1[Tidx] + (CoolData.reHeII1[Tidxp] - CoolData.reHeII1[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.reHeII1[i] = rate;
  
  rate = CoolData.reHeII2[Tidx] + (CoolData.reHeII2[Tidxp] - CoolData.reHeII2[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.reHeII2[i] = rate;
  
  rate = CoolData.reHeIII[Tidx] + (CoolData.reHeIII[Tidxp] - CoolData.reHeIII[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.reHeIII[i] = rate;
  
  rate = CoolData.brem[Tidx] + (CoolData.brem[Tidxp] - CoolData.brem[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.brem[i] = rate;
  
  rate = CoolData.hyd01k[Tidx] + (CoolData.hyd01k[Tidxp] - CoolData.hyd01k[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.hyd01k[i] = rate;
  
  rate = CoolData.h2k01[Tidx] + (CoolData.h2k01[Tidxp] - CoolData.h2k01[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.h2k01[i] = rate;
  
  rate = CoolData.vibh[Tidx] + (CoolData.vibh[Tidxp] - CoolData.vibh[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.vibh[i] = rate;
  
  rate = CoolData.roth[Tidx] + (CoolData.roth[Tidxp] - CoolData.roth[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.roth[i] = rate;
  
  rate = CoolData.rotl[Tidx] + (CoolData.rotl[Tidxp] - CoolData.rotl[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.rotl[i] = rate;
  
  rate = CoolData.GP99LowDensityLimit[Tidx] + (CoolData.GP99LowDensityLimit[Tidxp] - CoolData.GP99LowDensityLimit[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GP99LowDensityLimit[i] = rate;
  
  rate = CoolData.GP99HighDensityLimit[Tidx] + (CoolData.GP99HighDensityLimit[Tidxp] - CoolData.GP99HighDensityLimit[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GP99HighDensityLimit[i] = rate;
  
  rate = CoolData.HDlte[Tidx] + (CoolData.HDlte[Tidxp] - CoolData.HDlte[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.HDlte[i] = rate;
  
  rate = CoolData.HDlow[Tidx] + (CoolData.HDlow[Tidxp] - CoolData.HDlow[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.HDlow[i] = rate;
  
  rate = CoolData.HDcool[Tidx] + (CoolData.HDcool[Tidxp] - CoolData.HDcool[Tidx])*Tfac;
  for (i=0; i<5*CoolData.NumberOfTemperatureBins; i++)  CoolData.HDcool[i] = rate;
  
  rate = CoolData.cieco[Tidx] + (CoolData.cieco[Tidxp] - CoolData.cieco[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.cieco[i] = rate;
  
  rate = CoolData.GAHI[Tidx] + (CoolData.GAHI[Tidxp] - CoolData.GAHI[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GAHI[i] = rate;
  
  rate = CoolData.GAH2[Tidx] + (CoolData.GAH2[Tidxp] - CoolData.GAH2[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GAH2[i] = rate;
  
  rate = CoolData.GAHe[Tidx] + (CoolData.GAHe[Tidxp] - CoolData.GAHe[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GAHe[i] = rate;
  
  rate = CoolData.GAHp[Tidx] + (CoolData.GAHp[Tidxp] - CoolData.GAHp[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GAHp[i] = rate;
  
  rate = CoolData.GAel[Tidx] + (CoolData.GAel[Tidxp] - CoolData.GAel[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  CoolData.GAel[i] = rate;
  
  rate = RateData.k1[Tidx] + (RateData.k1[Tidxp] - RateData.k1[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k1[i] = rate;
  
  rate = RateData.k2[Tidx] + (RateData.k2[Tidxp] - RateData.k2[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k2[i] = rate;
  
  rate = RateData.k3[Tidx] + (RateData.k3[Tidxp] - RateData.k3[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k3[i] = rate;
  
  rate = RateData.k4[Tidx] + (RateData.k4[Tidxp] - RateData.k4[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k4[i] = rate;
  
  rate = RateData.k5[Tidx] + (RateData.k5[Tidxp] - RateData.k5[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k5[i] = rate;
  
  rate = RateData.k6[Tidx] + (RateData.k6[Tidxp] - RateData.k6[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k6[i] = rate;
  
  rate = RateData.k7[Tidx] + (RateData.k7[Tidxp] - RateData.k7[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k7[i] = rate;
  
  rate = RateData.k8[Tidx] + (RateData.k8[Tidxp] - RateData.k8[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k8[i] = rate;
  
  rate = RateData.k9[Tidx] + (RateData.k9[Tidxp] - RateData.k9[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k9[i] = rate;
  
  rate = RateData.k10[Tidx] + (RateData.k10[Tidxp] - RateData.k10[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k10[i] = rate;
  
  rate = RateData.k11[Tidx] + (RateData.k11[Tidxp] - RateData.k11[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k11[i] = rate;
  
  rate = RateData.k12[Tidx] + (RateData.k12[Tidxp] - RateData.k12[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k12[i] = rate;
  
  rate = RateData.k13[Tidx] + (RateData.k13[Tidxp] - RateData.k13[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k13[i] = rate;
  
  rate = RateData.k13dd[Tidx] + (RateData.k13dd[Tidxp] - RateData.k13dd[Tidx])*Tfac;
  for (i=0; i<7*CoolData.NumberOfTemperatureBins; i++)  RateData.k13dd[i] = rate;
  
  rate = RateData.k14[Tidx] + (RateData.k14[Tidxp] - RateData.k14[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k14[i] = rate;
  
  rate = RateData.k15[Tidx] + (RateData.k15[Tidxp] - RateData.k15[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k15[i] = rate;
  
  rate = RateData.k16[Tidx] + (RateData.k16[Tidxp] - RateData.k16[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k16[i] = rate;
  
  rate = RateData.k17[Tidx] + (RateData.k17[Tidxp] - RateData.k17[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k17[i] = rate;
  
  rate = RateData.k18[Tidx] + (RateData.k18[Tidxp] - RateData.k18[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k18[i] = rate;

  rate = RateData.k19[Tidx] + (RateData.k19[Tidxp] - RateData.k19[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k19[i] = rate;

  rate = RateData.k20[Tidx] + (RateData.k20[Tidxp] - RateData.k20[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k20[i] = rate;

  rate = RateData.k21[Tidx] + (RateData.k21[Tidxp] - RateData.k21[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k21[i] = rate;

  rate = RateData.k22[Tidx] + (RateData.k22[Tidxp] - RateData.k22[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k22[i] = rate;

  rate = RateData.k23[Tidx] + (RateData.k23[Tidxp] - RateData.k23[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k23[i] = rate;

  rate = RateData.k50[Tidx] + (RateData.k50[Tidxp] - RateData.k50[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k50[i] = rate;

  rate = RateData.k51[Tidx] + (RateData.k51[Tidxp] - RateData.k51[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k51[i] = rate;

  rate = RateData.k52[Tidx] + (RateData.k52[Tidxp] - RateData.k52[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k52[i] = rate;

  rate = RateData.k53[Tidx] + (RateData.k53[Tidxp] - RateData.k53[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k53[i] = rate;

  rate = RateData.k54[Tidx] + (RateData.k54[Tidxp] - RateData.k54[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k54[i] = rate;

  rate = RateData.k55[Tidx] + (RateData.k55[Tidxp] - RateData.k55[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k55[i] = rate;

  rate = RateData.k56[Tidx] + (RateData.k56[Tidxp] - RateData.k56[Tidx])*Tfac;
  for (i=0; i<CoolData.NumberOfTemperatureBins; i++)  RateData.k56[i] = rate;


  return SUCCESS;
}
