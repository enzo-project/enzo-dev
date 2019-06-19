/***********************************************************************
/
/  COMPUTE AND RETURN THE UNITS FOR CALCULATING HBAR
/
/  written by: Xinyu Li
/  date:       May, 2019
/  modified1:
/
/  PURPOSE:  Returns the cosmology units:
/
/         time:        utim = 1 / sqrt(4 * \pi * G * \rho_0
/         density:     urho = \rho_0 * (1+z)^3
/         length:      uxyz = (1 Mpc) * box / h
/         velocity:    uvel = uaye * uxyz / utim  (since u = a * dx/dt)
/    (*)  temperature: utem = m_H * \mu / k * uvel**2
/         a(t):        uaye = 1 / (1 + zri)
/
/           where:
/             box     - size of simulation box in Mpc/h
/             zri     - initial redshift (start of simulation)
/             \rho_0  = 3*\Omega_0*H_0^2/(8*\pi*G)
/             Omega_0 - the fraction of non-relativistic matter at z=0
/
/   NOTE: the units are calculated at a=1 (initial redshift)
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"
#include "units.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "TopGridData.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int QuantumCosmologyGetUnits(float *DensityUnits, float *LengthUnits,
          float *TemperatureUnits, float *TimeUnits,
          float *VelocityUnits, FLOAT Time)
{
 
  /* From the time, compute the current redshift. */
 
  FLOAT a, dadt;
  if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    ENZO_FAIL("Error in ComputeExpansionFactor.\n");

  }
 
  /* Compute the current redshift (remember a(init) = 1). */
 
  FLOAT CurrentRedshift = (1 + InitialRedshift)/a - 1;
 
  /* Determine the units. */
 
  *DensityUnits     = 1.8788e-29*OmegaMatterNow*POW(HubbleConstantNow,2)*
                      POW(1 + CurrentRedshift,3);
 
  *LengthUnits      = 3.085678e24*ComovingBoxSize/HubbleConstantNow/
                      (1 + InitialRedshift);
 
  *TemperatureUnits = 1.81723e6*POW(ComovingBoxSize,2)*OmegaMatterNow*
                      (1 + InitialRedshift);
 
  *TimeUnits        = 2.519445e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                      POW(1 + InitialRedshift,FLOAT(1.5));
 
  *VelocityUnits    = 1.22475e7*ComovingBoxSize*sqrt(OmegaMatterNow)*
                      sqrt(1 + InitialRedshift);
 
  return SUCCESS;
}

int QuantumGetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time)
{
  /* If using cosmology, get cosmology units */
  if (ComovingCoordinates) {
    if (QuantumCosmologyGetUnits(DensityUnits, LengthUnits, TemperatureUnits,
                          TimeUnits, VelocityUnits, Time) == FAIL) {
      ENZO_FAIL("Error in QuantumCosmologyGetUnits.\n");
    }
  }
  else {
      /* Determine the units. */
      *DensityUnits = GlobalDensityUnits;
      *MassUnits = GlobalMassUnits;
      *LengthUnits      = GlobalLengthUnits;
      *TemperatureUnits = mh*pow(GlobalLengthUnits/GlobalTimeUnits,2)/kboltz;   //K
      *TimeUnits        = GlobalTimeUnits;
      *VelocityUnits    = GlobalLengthUnits/GlobalTimeUnits; //cms-1
  
    }
  return SUCCESS;
}

int QuantumGetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, FLOAT Time)

{
  /* If using cosmology, get cosmology units */
  if (ComovingCoordinates) {
    //    fprintf(stderr, "Using CosmologyCoordinates.\n");
    if (QuantumCosmologyGetUnits(DensityUnits, LengthUnits, TemperatureUnits,

                          TimeUnits, VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in QuantumCosmologyGetUnits.\n");
      return FAIL;
    }
  }
  else {
      /* Determine the units. */
      *DensityUnits = GlobalDensityUnits;
      *LengthUnits      = GlobalLengthUnits;
      *TemperatureUnits = mh*pow(GlobalLengthUnits/GlobalTimeUnits,2)/kboltz;   //K
      *TimeUnits        = GlobalTimeUnits;
      *VelocityUnits    = GlobalLengthUnits/GlobalTimeUnits; //cms-1
    }
  return SUCCESS;
}




