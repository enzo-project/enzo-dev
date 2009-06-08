#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"

/* These are the different types of baryon fields. */

enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
                 FieldUndefined};

#define FieldTypeIsDensity(A) (((A) >= TotalEnergy && (A) <= Velocity3) ? FALSE : TRUE)

/* These are the different types of fluid boundary conditions. */

enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated, 
				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro};

/* Define a float/int union. */

union float_int {
  int ival;
  float fval;
  FLOAT FVAL;
};

#endif
