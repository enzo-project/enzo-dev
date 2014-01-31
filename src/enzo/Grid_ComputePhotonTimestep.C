/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    dt   - timestep
/
************************************************************************/

// Compute the timestep from all the constrains for this grid.
//
// Somebody fix the error handling in this routine! please.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferParameters.h"

/* function prototypes */

int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
extern "C" void FORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
			     hydro_method *ihydro, float *C2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
                             float *vgz, float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);


float grid::ComputePhotonTimestep()
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;

  this->DebugCheck((char*) "ComputeTimeStep");

  /* initialize */

  float dt, dtTemp;
  float dtBaryons      = huge_number;
  float dtViscous      = huge_number;
  float dtParticles    = huge_number;
  float dtExpansion    = huge_number;
  float dtAcceleration = huge_number;
  int dim, i, result;

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, TimeUnits, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Compute the field size. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  float afloat = float(a);

  /* 1) Compute Courant condition for baryons. */

  if (NumberOfBaryonFields > 0) {

    /* Find fields: density, total energy, velocity1-3. */
  
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }

    /* Compute the pressure. */

    float *pressure_field = new float[size];
    result = this->ComputePressure(Time, pressure_field);

    if (result == FAIL) {
      fprintf(stderr, "Error in grid->ComputePressure.\n");
      exit(EXIT_FAILURE);
    }

#ifdef UNUSED
    int Zero[3] = {0,0,0}, TempInt[3] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */

    /* Call fortran routine to do calculation. */
 
    FORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
                               GridDimension+2,
//                        Zero, TempInt, Zero+1, TempInt+1, Zero+2, TempInt+2,
                          GridStartIndex, GridEndIndex,
                               GridStartIndex+1, GridEndIndex+1,
                               GridStartIndex+2, GridEndIndex+2,
			       &HydroMethod, &ZEUSQuadraticArtificialViscosity,
                          CellWidth[0], CellWidth[1], CellWidth[2],
                               GridVelocity, GridVelocity+1, GridVelocity+2,
                               &Gamma, &PressureFree, &afloat,
                          BaryonField[DensNum], pressure_field,
                               BaryonField[Vel1Num], BaryonField[Vel2Num],
                               BaryonField[Vel3Num], &dtBaryons, &dtViscous);

    /* Clean up */

    delete [] pressure_field;

    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */

    dtBaryons *= CourantSafetyNumber;
    
  }

  /* 2) Calculate dt from particles. */

  if (NumberOfParticles > 0) {
    
    /* Compute dt constraint from particle velocities. */

    for (dim = 0; dim < GridRank; dim++) {
      float dCell = CellWidth[dim][0]*a;
      for (i = 0; i < NumberOfParticles; i++) {
        dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
	dtParticles = min(dtParticles, dtTemp);
      }
    }

    /* Multiply resulting dt by ParticleCourantSafetyNumber. */

    dtParticles *= ParticleCourantSafetyNumber;

  }

  /* 3) Find dt from expansion. */

  if (ComovingCoordinates)
    if (CosmologyComputeExpansionTimestep(Time, &dtExpansion) == FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(FAIL);
    }

  /* 4) Calculate minimum dt due to acceleration field (if present). */

  if (SelfGravity) {
    for (dim = 0; dim < GridRank; dim++)
      if (AccelerationField[dim])
	for (i = 0; i < size; i++) {
	  dtTemp = sqrt(CellWidth[dim][0]/
			fabs(AccelerationField[dim][i])+tiny_number);
	  dtAcceleration = min(dtAcceleration, dtTemp);
	}
    if (dtAcceleration != huge_number)
      dtAcceleration *= 0.5;
  }

  /* 5) calculate minimum timestep */

  dt = min(dtBaryons, dtParticles);
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);

  /* 6) If star formation (Pop III for now), set a minimum timestep */

#ifdef UNUSED
  float mindtNOstars;  // Myr
  const int NumberOfStepsInLifetime = 5;
  float dtStar = huge_number;

  if (STARFEED_METHOD(POP3_STAR))
    mindtNOstars = 3;  // Myr
  if (STARFEED_METHOD(STAR_CLUSTER))
    mindtNOstars = 10;  // Myr

  if (STARFEED_METHOD(POP3_STAR) && STARFEED_METHOD(STAR_CLUSTER))
    if (G_TotalNumberOfStars > 0 && minStarLifetime < 1e6)
      dtStar = minStarLifetime/NumberOfStepsInLifetime;
    else
      dtStar = 3.1557e13*mindtNOstars/TimeUnits;

  dt = min(dt, dtStar);
#endif
  
  /* 7) If using radiation pressure, calculate minimum dt */

  float dtRadPressure = huge_number;
  float absVel, absAccel;

  if (RadiationPressure && RadiativeTransfer) {

    int RPresNum1, RPresNum2, RPresNum3;
    if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) 
	== FAIL) {
      ENZO_FAIL("Error in IdentifyRadiationPressureFields.\n");
    }

    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++) {
	dtTemp = sqrt(CellWidth[dim][0] / (fabs(BaryonField[RPresNum1+dim][i])+
					   tiny_number));
	dtRadPressure = min(dtRadPressure, dtTemp);
      }
    
    if (dtRadPressure < huge_number)
      dtRadPressure *= 0.5;

    dt = min(dt, dtRadPressure);

  } /* ENDIF RadiationPressure */

  /* 8) Safety Velocity to limit timesteps */

  float dtSafetyVelocity = huge_number;
  if (TimestepSafetyVelocity > 0)
    dtSafetyVelocity = a*CellWidth[0][0] / 
      (TimestepSafetyVelocity*1e5 / VelocityUnits);    // parameter in km/s

  dt = min(dt, dtSafetyVelocity);

  /* 9) If we're calculating the timestep for RT, limit it (doesn't
     affect hydro dt). */

  // parameter in km/s
  float dx_level, dx_ratio;
  float dtPhotonSafety = tiny_number;

  if (RadiativeTransferTimestepVelocityLimit > 0)
    dtPhotonSafety = a*CellWidth[0][0] / 
      (RadiativeTransferTimestepVelocityLimit*1e5 / VelocityUnits);
  if (RadiativeTransferTimestepVelocityLevel >= 0) {
    dx_level = TopGridDx[0] * POW(RefineBy, -RadiativeTransferTimestepVelocityLevel);
    dx_ratio = dx_level / CellWidth[0][0];
    if (dx_ratio > 1)
      dtPhotonSafety *= dx_ratio;
  }
  dt = max(dt, CourantSafetyNumber*dtPhotonSafety);

  /* Debugging info. */

//  if (debug || NumberOfProcessors > 1) {
//  if (debug) {
//    printf("ComputeTimeStep = %"FSYM" (", dt);
//    if (NumberOfBaryonFields > 0)
//      printf("Bar = %"GSYM" ", dtBaryons);
//    if (HydroMethod == Zeus_Hydro)
//      printf("Vis = %"FSYM" ", dtViscous);
//    if (ComovingCoordinates)
//      printf("Exp = %"FSYM" ", dtExpansion);
//    if (dtAcceleration != huge_number)
//      printf("Acc = %"FSYM" ", dtAcceleration);
//    if (NumberOfParticles)
//      printf("Part = %"FSYM" ", dtParticles);
//#ifdef TRANSFER
//    if (RadiationPressure && RadiativeTransfer && dtRadPressure < 100)
//      printf("Rad = %"GSYM" ", dtRadPressure);
//    if (dtSafetyVelocity != huge_number)
//      printf("Saf = %"GSYM" ", dtSafetyVelocity); 
//#endif /* TRANSFER */
//    printf(")\n");
//  }

  return dt;
}
