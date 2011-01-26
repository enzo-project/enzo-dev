/***********************************************************************
/
/  GRID CLASS (ADD FIXED ACCELERATION TO ACCEL FIELDS)
/
/  written by: Greg Bryan
/  date:       June, 1996
/  modified1:  Brian O'Shea  (fixed NFW profile stuff)
/  date:       July  2009
/
/  PURPOSE: Certain problems required external acceleration fields.
/    This routine adds them to the existing self-gravitating fields, or
/     creates new fields if necessary.
/    There is currently support for:
/      1) A uniform gravitational field in one of the orthogonal directions
/         for grid cells only
/      2) A 3D point source field for grid cells only
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"
 
/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
       	      float *TemperatureUnits, float *TimeUnits,
       	      float *VelocityUnits, double *MassUnits, FLOAT Time);
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
int grid::ComputeAccelerationFieldExternal()
{
 
  /* Return if this does not concern us */
  if (!(UniformGravity || PointSourceGravity || ExternalGravity)) return SUCCESS;

  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  int dim, i, j, k, size = 1;
 
  LCAPERF_START("grid_ComputeAccelerationFieldExternal");

  /* Compute field size (in floats). */
 
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Check if acceleration field exists.  If not create it and zero it. */
 
  if (AccelerationField[0] == NULL)
    for (dim = 0; dim < GridRank; dim++) {
      AccelerationField[dim] = new float[size];
      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = 0;
    }
 

  /* BWO:  if we're using the NFW halo with ProblemType = 31 ("GalaxySimulation"),
     the PointSourceGravitationalConstant is assumed to be the mass in 
     actually a mass in code units (that needs to be converted to CGS).  If this
     problem type is used, make this conversion.  Otherwise, don't worry about it. */
  double MassUnitsDouble;
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1, AccelUnits = 1;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  AccelUnits = LengthUnits/TimeUnits/TimeUnits;
  if(PointSourceGravity == 2 && ProblemType == 31){
    if(ComovingCoordinates){
      MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits),3.0);
    } else {
      MassUnitsDouble = double(MassUnits);
    }
  } else {
    MassUnitsDouble = 1.0;
  }


  /* -----------------------------------------------------------------
     Point Source gravity
     ----------------------------------------------------------------- */
 
  if (PointSourceGravity) {
 
    FLOAT a = 1.0, accel, dadt, radius, rcubed, rsquared, 
      xpos, ypos = 0.0, zpos = 0.0, rcore,x ;
 
    /* Compute adot/a at time = t+1/2dt (time-centered). */
 
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	  == FAIL) {
		ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
      }

    /* Loop over grid, adding acceleration to field. */

    if (PointSourceGravity == 1 || PointSourceGravity == 2)
      for (dim = 0; dim < GridRank; dim++) {
	int n = 0;
 
	for (k = 0; k < GridDimension[2]; k++) {
	  if (GridRank > 2)
	    zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] -
	      PointSourceGravityPosition[2];
	  if (dim == 2 && HydroMethod == Zeus_Hydro)
	    zpos -= 0.5*CellWidth[2][k];
	  
	  for (j = 0; j < GridDimension[1]; j++) {
	    if (GridRank > 1)
	      ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] -
		PointSourceGravityPosition[1];
	    if (dim == 1 && HydroMethod == Zeus_Hydro)
	      ypos -= 0.5*CellWidth[1][j];
	    
	    for (i = 0; i < GridDimension[0]; i++, n++) {
	      xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] -
		PointSourceGravityPosition[0];
	      if (dim == 0 && HydroMethod == Zeus_Hydro)
		xpos -= 0.5*CellWidth[0][i];

	      /* Compute distance from center. */
 
	      rsquared = xpos*xpos + ypos*ypos + zpos*zpos;
	      rcubed = POW(rsquared, 1.5);
 
	      if (PointSourceGravity == 1) {
		
		/* (1) Point Source:
		   Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm(). 
		   (i.e. 1/a^2 * a = 1/a). */
		
		rcore = max(0.1*CellWidth[0][0], PointSourceGravityCoreRadius);
		accel = min(PointSourceGravityConstant/((rsquared)*POW(rsquared,0.5)*a),
			    PointSourceGravityConstant/(rcore*rcore*POW(rsquared,0.5)*a));
		
	      } else if (PointSourceGravity == 2) {
		  
		/* (2) NFW Profile: CoreRadius and Constant are both in code units 
		   if ProblemType == 31 (Galaxy simulation), otherwise they're in CGS.
		   Need to convert the core radius to code units and the gravity constant to
		   CGS.  (BWO, July 2009)  */

		radius = sqrt(rsquared);
		if(ProblemType == 31){
		  rcore = PointSourceGravityCoreRadius;  // already in code units
		} else {
		  rcore = PointSourceGravityCoreRadius/LengthUnits;  // convert from CGS to code
		}

		FLOAT x = radius/rcore;

		// BWO, July 2009: MassUnitsDouble is CGS mass units if ProblemType == 31,
		// and 1.0 otherwise.
		accel = GravConst*PointSourceGravityConstant*MassUnitsDouble*
		  ((log(1.0+x  )-x  /(1.0+x  )) /
		   (log(1.0+1.0)-1.0/(1.0+1.0))) / 
		  POW(radius*LengthUnits, 2.0) / AccelUnits;

		accel = accel/radius;  // this radius normalizes the multiplication by 
		                         // xpos,ypos,zpos done below

	      } else {
		/* this is only reached if there are two types of point sources - 
		   when you add a new one, this changes */
		ENZO_FAIL("should never get here! in Grid::ComputeAccelFieldExternal");
	      }

	      /* Apply force. */
	
	      if (dim == 0)
		AccelerationField[0][n] -= accel*xpos;
	      if (dim == 1)
		AccelerationField[1][n] -= accel*ypos;
	      if (dim == 2)
		AccelerationField[2][n] -= accel*zpos;

	    } // end: loop over i
	  } // end: loop over j
	} // end: loop over k
      } // end: loop over dims
   
    /* DO PARTICLES HERE! */

    if (NumberOfParticles > 0 && GridRank != 3) {
        ENZO_FAIL("PointSourceGravity assumes 3D");
    }
      
    if (PointSourceGravity == 1 || PointSourceGravity == 2)
      for (i = 0; i < NumberOfParticles; i++) {
	
	/* Compute vector between particle (advanced by 1/2 step) and
	   gravity center. */
      
	xpos = ParticlePosition[0][i] + 0.5*dtFixed*ParticleVelocity[0][i]/a - 
	  PointSourceGravityPosition[0];
	ypos = ParticlePosition[1][i] + 0.5*dtFixed*ParticleVelocity[1][i]/a - 
	  PointSourceGravityPosition[1];
	zpos = ParticlePosition[2][i] + 0.5*dtFixed*ParticleVelocity[2][i]/a - 
	  PointSourceGravityPosition[2];
	    
	/* Compute distance from center. */

	rsquared = xpos*xpos + ypos*ypos + zpos*zpos;
      
	/* (1) is a real (softened) point-source, (2) is NFW profile */
	
	if (PointSourceGravity == 1) {

	  /* (1) Point Source:
	     Multiply by a(t) to offset the 1/a(t) in ComovingAccelTerm(). 
	     (i.e. 1/a^2 * a = 1/a). */

	  rcore = max(0.1*CellWidth[0][0], PointSourceGravityCoreRadius);
	  accel = min(PointSourceGravityConstant/(rsquared*POW(rsquared,0.5)*a),
		      PointSourceGravityConstant/(rcore*rcore*POW(rsquared,0.5)*a));

	}  // if (PointSourceGravity == 1)

	if (PointSourceGravity == 2) {

	  /* (2) NFW Profile: assume CoreRadius is rs in cm and Constant
	     is mass within rs in g. */

	  
	  radius = sqrt(rsquared);
	  rcore = PointSourceGravityCoreRadius/LengthUnits;
	  FLOAT x = radius/rcore;
	  accel = GravConst*PointSourceGravityConstant*
	    ((log(1.0+x  )-x  /(1.0+x  )) /
	     (log(1.0+1.0)-1.0/(1.0+1.0))) / 
	    POW(radius*LengthUnits, 2) / AccelUnits;
	  accel = accel/radius;  // this radius normalizes the multiplication by xpos,ypos,zpos done below


	} // if (PointSourceGravity == 2)
	

	/* Apply force. */
	
	ParticleAcceleration[0][i] -= accel*xpos;
	ParticleAcceleration[1][i] -= accel*ypos;
	ParticleAcceleration[2][i] -= accel*zpos;

      } // end: loop over number of particles
      
  } // end: if (PointSourceGravity)


  /* -----------------------------------------------------------------
     ExternalGravity 1: another similar way for a NFW profile
     ----------------------------------------------------------------- */

  if (ExternalGravity == 1) {
    
    /* Specify NFW parameters by hand here 
       Should move to parameter file in the future */
    
    double rhoc = HaloCentralDensity,
      c = HaloConcentration, 
      rvir = HaloVirialRadius;
    FLOAT xc = 0.5, yc = 0.5, zc = 0.5;

    double rs = rvir / c;
    double Mvir = 4.0*M_PI*rhoc*pow(rs,3)*(log(1.0+c)-c/(1.0+c));
    
    float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
      TimeUnits = 1.0, VelocityUnits = 1.0;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	     &TimeUnits, &VelocityUnits, &MassUnits, Time);
    double AccelerationUnits = LengthUnits / pow(TimeUnits,2);
    double CGSGravConst = 6.672e-8;

    printf("rhoc=%g, rvir=%g, Mvir=%g\n", rhoc, rvir, Mvir/1.989e33);
    
    FLOAT x, y, z, xpos, ypos, zpos, r;
    int n = 0;
    double x1, M, g;
    for (dim = 0; dim < GridRank; dim++) {
      n = 0;
    for (int k = 0; k < GridDimension[2]; k++) {
      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
      zpos = z - zc;

      if (dim == 2 && HydroMethod == Zeus_Hydro) 
	zpos -= 0.5*CellWidth[2][k];

      for (int j = 0; j < GridDimension[1]; j++) {
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	ypos = y - yc;

	if (dim == 1 && HydroMethod == Zeus_Hydro) 
	  ypos -= 0.5*CellWidth[1][j];

	for (int i = 0; i < GridDimension[0]; i++, n++) {
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  xpos = x - xc;

	  if (dim == 0 && HydroMethod == Zeus_Hydro) 
	    xpos -= 0.5*CellWidth[0][i];

	  r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
	  r = max(r, CellWidth[0][0]);
	  
	  if (r < rvir/LengthUnits) {
	    x1 = r*LengthUnits/rs;
	    M = 4.0*M_PI*rhoc*pow(rs,3)*(log(1.0+x1)-x1/(1.0+x1));
	  }
	  else {
	    M = Mvir;
	  }
	  g = CGSGravConst*M/pow(r*LengthUnits,2);
	  g /= AccelerationUnits;
	  if (dim == 0) { 
	    AccelerationField[0][n] += -g*xpos/r;
	  }
	  if (dim == 1) {
	    AccelerationField[1][n] += -g*ypos/r;
	  }
	  if (dim == 2) { 
	    AccelerationField[2][n] += -g*zpos/r;
	  }
	}
      }
    }
    }

    for (int i = 0; i < NumberOfParticles; i++) {
      x = ParticlePosition[0][i];
      y = ParticlePosition[1][i];
      z = ParticlePosition[2][i];
      xpos = x - xc;
      ypos = y - yc;
      zpos = z - zc;
      r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
      r = max(r, CellWidth[0][0]);

      if (r < rvir/LengthUnits) {
	x1 = r*LengthUnits/rs;
	M = 4.0*M_PI*rhoc*pow(rs,3)*(log(1.0+x1)-x1/(1.0+x1));
      }
      else {
	M = Mvir;
      }
      g = CGSGravConst*M/pow(r*LengthUnits,2);
      g /= AccelerationUnits;

      ParticleAcceleration[0][i] += -g*xpos/r;
      ParticleAcceleration[1][i] += -g*ypos/r;
      ParticleAcceleration[2][i] += -g*zpos/r;
    }
    

  }


  /* -----------------------------------------------------------------
     Uniform gravity field
     ----------------------------------------------------------------- */
 
  if (UniformGravity) {
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* Set constant for this dimension. */
 
      float Constant = 0.0;
      if (dim == UniformGravityDirection)
	Constant = UniformGravityConstant;
	
      /* Set field. */
 
      for (i = 0; i < size; i++)
	AccelerationField[dim][i] = Constant;
 
    } // loop over dims
 
  } // end: if (UniformGravity)
 
  LCAPERF_STOP("grid_ComputeAccelerationFieldExternal");
  return SUCCESS;
}
 
