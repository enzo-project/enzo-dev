/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A ZELDOVICH PANCAKE)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
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
#include "CosmologyParameters.h"
 
#define TOLERANCE 1.0e-6
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::ZeldovichPancakeInitializeGrid(int  ZeldovichPancakeDirection,
				      float ZeldovichPancakeCentralOffset,
				      float ZeldovichPancakeOmegaBaryonNow,
				      float ZeldovichPancakeOmegaCDMNow,
				      float ZeldovichPancakeCollapseRedshift,
				      float ZeldovichPancakeInitialTemperature,
				      float ZeldovichPancakeInitialUniformBField[]
					 )
{
  /* declarations */
 
  float Amplitude, AmplitudeVel, kx, xLagrange, xEulerian, xEulerianOld;
  int   dim, field, i, index;
  const float Pi = 3.14159;
 
  /* error check */
 
  if (ZeldovichPancakeDirection < 0 || ZeldovichPancakeDirection >= GridRank) {
    ENZO_FAIL("ZeldovichPancakeDirection is improperly set.\n");
  }
  if (ZeldovichPancakeOmegaCDMNow != 0) {
    ENZO_FAIL("Dark matter not yet supported.\n");
  }
 
  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (MaxVelocityIndex>1){
    FieldType[NumberOfBaryonFields++] = Velocity2;
  }
  if (MaxVelocityIndex > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  int iTE = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Get the cosmology units so we can convert temperature later. */
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits,
	       InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  /* Determine the size of the fields. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate space for the fields. */
 
  this->AllocateGrids();
 
  /* Find the stride between zones along the pancake direction. */
 
  int Divisor = 1;
  for (dim = 0; dim < ZeldovichPancakeDirection; dim++)
    Divisor *= GridDimension[dim];
 
  int NumberOfZones = GridEndIndex[ZeldovichPancakeDirection] -
                      GridStartIndex[ZeldovichPancakeDirection] + 1;
 
  /* Compute the amplitude of perturbations. */
 
  Amplitude    = (1+ZeldovichPancakeCollapseRedshift) / (1+InitialRedshift);
  AmplitudeVel = -sqrt(2.0/3.0)*(1.0+ZeldovichPancakeCollapseRedshift) /
                 ((1.0+InitialRedshift) * 2.0*Pi);
  kx           = 2*Pi/float(NumberOfZones);
 
  /* set density, total energy and velocity in problem dimension */
 
  for (i = 0; i < size; i++) {
 
    /* Determine the index along the pancake direction and convert this
       into a number from -0.5NumberOfZones to +0.5NumberOfZones. */
 
    index     = (i/Divisor % GridDimension[ZeldovichPancakeDirection]) -
                GridStartIndex[ZeldovichPancakeDirection];
    index     = (index+NumberOfZones) % NumberOfZones;
    xLagrange = float(index) + (0.5-ZeldovichPancakeCentralOffset) -
                0.5*float(NumberOfZones);
 
    /* Convert this Lagrangean position into an Eulerian one using
       a Newton-Raphson method. */
 
    xEulerian    = xLagrange;
    xEulerianOld = FLOAT_UNDEFINED;
    while (fabs((xEulerian-xEulerianOld)/xEulerian) > TOLERANCE) {
      xEulerianOld = xEulerian;
      xEulerian += (xLagrange - xEulerian + Amplitude*sin(kx*xEulerian)/kx)
	          /(1                     - Amplitude*cos(kx*xEulerian)   );
    }
 
    /* Set density. */
    // correct Zeldovich test: 
    BaryonField[0][i] = ZeldovichPancakeOmegaBaryonNow/
	  (1 - Amplitude*cos(kx*xEulerian));
    // terribly fudge since the folks that do B field tests
    // did not set up the density fields consistently ...
    // MATCH RYU 1993 et al paper ..... 
    //      BaryonField[0][i] = 1.;

    /* Set total energy, gas energy. */
 
    BaryonField[iTE][i] = ZeldovichPancakeInitialTemperature/TemperatureUnits *
          POW(BaryonField[0][i]/ZeldovichPancakeOmegaBaryonNow,Gamma-1)
	                  / (Gamma-1);
    if (DualEnergyFormalism)
      BaryonField[iTE+1][i] = BaryonField[iTE][i];
    if (HydroMethod != Zeus_Hydro)
      BaryonField[iTE][i]  += 0.5*POW(AmplitudeVel*sin(kx*xEulerian),2);
 
    /* Set velocities (some may be set to zero later -- see below). */
 
    if (HydroMethod == Zeus_Hydro) xEulerian -= 0.5; // only approximate
    BaryonField[vel][i]     = AmplitudeVel * sin(kx*xEulerian);
    if (GridRank > 1 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
      BaryonField[vel+1][i] = AmplitudeVel * sin(kx*xEulerian);
    if (GridRank > 2 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
      BaryonField[vel+2][i] = AmplitudeVel * sin(kx*xEulerian);

    if ( UseMHD ) {
      BaryonField[iBx  ][i] = ZeldovichPancakeInitialUniformBField[0];
      BaryonField[iBy  ][i] = ZeldovichPancakeInitialUniformBField[1];
      BaryonField[iBz  ][i] = ZeldovichPancakeInitialUniformBField[2];
      BaryonField[ietot][i] += 0.5*(BaryonField[iBx][i] * BaryonField[iBx][i]+
				    BaryonField[iBy][i] * BaryonField[iBy][i]+
				    BaryonField[iBz][i] * BaryonField[iBz][i])/
	BaryonField[iden][i];
    }
    if ( UseMHDCT ){
        MagneticField[0][i] = ZeldovichPancakeInitialUniformBField[0];
        MagneticField[1][i] = ZeldovichPancakeInitialUniformBField[1];
        MagneticField[2][i] = ZeldovichPancakeInitialUniformBField[2];
    }
    if (HydroMethod == MHD_RK) {
      BaryonField[iPhi ][i] = 0.0;
    }
  }
 

  /* set transverse velocities (i.e. erase any incorrectly set velocities). */
 
  for (field = vel; field < vel+((HydroMethod == MHD_RK  || (HydroMethod == HD_RK) ) ? 3 : GridRank); field++)
    if (field != ZeldovichPancakeDirection+vel)
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0.0;
 
  return SUCCESS;
}
