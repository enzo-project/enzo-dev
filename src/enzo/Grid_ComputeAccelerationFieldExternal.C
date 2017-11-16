/***********************************************************************
/
/  GRID CLASS (ADD FIXED ACCELERATION TO ACCEL FIELDS)
/
/  written by: Greg Bryan
/  date:       June, 1996
/  modified1:  Brian O'Shea  (fixed NFW profile stuff)
/  date:       July  2009
/  modified2:  Elizabeth Tasker (added in call for external potential)
/  date:      September 2011
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
  if (!(UniformGravity || PointSourceGravity || DiskGravity || ExternalGravity)) return SUCCESS;

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
    if( SelfGravity == 0 ){
      for (dim = 0; dim < GridRank; dim++) {
        for (i = 0; i < size; i++){
          AccelerationField[dim][i] = 0;
        }
      }
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
 
  if (PointSourceGravity > 0) {
 
    FLOAT a = 1.0, accel, dadt, radius, rcubed, rsquared, 
      xpos, ypos = 0.0, zpos = 0.0, rcore,x ;
 
    /* Compute adot/a at time = t+1/2dt (time-centered). */
 
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt)
	  == FAIL) {
		ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
      }

    /* Loop over grid, adding acceleration to field. */

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
              //Yuan, Aug 2011: Add BCG and SMBH potential if ProblemType == 108
              if(ProblemType == 108){
                accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                ClusterSMBHBCG*POW((POW(POW(radius*LengthUnits/(1.0e-3*Mpc),0.5975)/3.206e-7,0.9) +
                POW(POW(radius*LengthUnits/(1.0e-3*Mpc), 1.849)/1.861e-6, 0.9)), -1.0/0.9) / AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass / POW(radius*LengthUnits, 2) / AccelUnits; // + BCG + BH mass
                  /*Bondi*/
                if(ClusterSMBHCalculateGasMass == 4){
                 accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                ClusterSMBHBCG*POW((POW(POW(radius*LengthUnits/(1.0e-3*Mpc),0.5975)/3.206e-7,0.9) +
                POW(POW(radius*LengthUnits/(1.0e-3*Mpc), 1.849)/1.861e-6, 0.9)), -1.0/0.9) / AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass/POW(radius*LengthUnits - 2.0*GravConst*SolarMass/POW(clight,2), 2)/ AccelUnits;
                 }
                  /*Elliptical Galaxy Fixed Gravity*/
                if(EllipticalGalaxyRe > 0.001){  
                 accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                GravConst*(ClusterSMBHBCG*SolarMass*1.0e11)/POW(radius*LengthUnits+EllipticalGalaxyRe*1.0e-3*Mpc/1.8153, 2)/AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass/POW(radius*LengthUnits - 2.0*GravConst*SolarMass/POW(clight,2), 2)/ AccelUnits;
                 }
              }
              accel = accel/radius;  // this radius normalizes the multiplication by
	      // xpos,ypos,zpos done below

	    }  else if (PointSourceGravity == 3) {
	      // (3) Isothermal sphere along a line at the cylindrical radius, a,  
	      //    
	      FLOAT a2 = PointSourceGravityCoreRadius*PointSourceGravityCoreRadius;
	      /*	      if (GridRank > 1)
		a2 = max(a2, ypos*ypos);
	      if (GridRank > 2)
		a2 += zpos*zpos;
	      */
	      //	      accel = PointSourceGravityConstant/(xpos*xpos + a2);
	      accel = PointSourceGravityConstant/(rsquared + a2);
	      //	      fprintf(stderr, "%g %g %g\n", xpos, accel, AccelUnits);
	    }  else {
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
      
    if (PointSourceGravity > 0)
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
	  
	  if(ProblemType == 31){
	    rcore = PointSourceGravityCoreRadius;  // already in code units                                                      
	  } else {
	    rcore = PointSourceGravityCoreRadius/LengthUnits;  // convert from CGS to code                                       
	  }
	  FLOAT x = radius/rcore;
	  accel = GravConst*PointSourceGravityConstant*MassUnitsDouble*
	    ((log(1.0+x  )-x  /(1.0+x  )) /
	     (log(1.0+1.0)-1.0/(1.0+1.0))) / 
	    POW(radius*LengthUnits, 2) / AccelUnits;
              /*Yuan, Aug 2011: Add BCG and SMBH potential if ProblemType == 108*/
              if(ProblemType == 108){
                accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                ClusterSMBHBCG*POW((POW(POW(radius*LengthUnits/(1.0e-3*Mpc),0.5975)/3.206e-7,0.9) +
                POW(POW(radius*LengthUnits/(1.0e-3*Mpc), 1.849)/1.861e-6, 0.9)), -1.0/0.9) / AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass / POW(radius*LengthUnits, 2) / AccelUnits; // + BCG + BH mass
                  /*Bondi*/
                if(ClusterSMBHCalculateGasMass == 4){
                 accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                ClusterSMBHBCG*POW((POW(POW(radius*LengthUnits/(1.0e-3*Mpc),0.5975)/3.206e-7,0.9) +
                POW(POW(radius*LengthUnits/(1.0e-3*Mpc), 1.849)/1.861e-6, 0.9)), -1.0/0.9) / AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass/POW(radius*LengthUnits - 2.0*GravConst*SolarMass/POW(clight,2), 2)/ AccelUnits;
                 }
                  /*Elliptical Galaxy Fixed Gravity*/
                if(EllipticalGalaxyRe > 0.001){
                 accel = GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x  )-x  /(1.0+x  )) /
                 (log(1.0+1.0)-1.0/(1.0+1.0))) /
                POW(radius*LengthUnits, 2.0) / AccelUnits +
                GravConst*(ClusterSMBHBCG*SolarMass*1.0e11)/POW(radius*LengthUnits+EllipticalGalaxyRe*1.0e-3*Mpc/1.8153, 2)/AccelUnits +
                GravConst*SolarMass*ClusterSMBHMass/POW(radius*LengthUnits - 2.0*GravConst*SolarMass/POW(clight,2), 2)/ AccelUnits;
                 }
              }
	  accel = accel/radius;  // this radius normalizes the multiplication by xpos,ypos,zpos done below


	} // if (PointSourceGravity == 2)
	

	/* Apply force. */
	
	ParticleAcceleration[0][i] -= accel*xpos;
	ParticleAcceleration[1][i] -= accel*ypos;
	ParticleAcceleration[2][i] -= accel*zpos;

      } // end: loop over number of particles
      
  } // end: if (PointSourceGravity)


/*-----------------------------------------------------------------------
 *     DiskGravity
 *     Reference: Burkert 1995, Mori & Burkert 2000
 *------------------------------------------------------------------------*/

  if (DiskGravity > 0) {

    double accel, radius, rsquared, xpos, ypos = 0, zpos = 0, rcore,rcyl;
    FLOAT dadt, a = 1;
    float AngularMomentumx, AngularMomentumy, AngularMomentumz;
    float MSDisk, SDiskScaleHeightR, SDiskScaleHeightz, MBulge, rBulge,
      rDMConst, densDMConst;

    AngularMomentumx = DiskGravityAngularMomentum[0];
    AngularMomentumy = DiskGravityAngularMomentum[1];
    AngularMomentumz = DiskGravityAngularMomentum[2];
    MSDisk = DiskGravityStellarDiskMass;
    SDiskScaleHeightR = DiskGravityStellarDiskScaleHeightR;
    SDiskScaleHeightz = DiskGravityStellarDiskScaleHeightz;
    MBulge = DiskGravityStellarBulgeMass;
    rBulge = DiskGravityStellarBulgeR;
    rDMConst = DiskGravityDarkMatterR;
    densDMConst = DiskGravityDarkMatterDensity;

    /* Compute adot/a at time = t+1/2dt (time-centered). */
    float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
          VelocityUnits=1, AccelUnits=1;
    double MassUnits=1;

    if (ComovingCoordinates) {
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) == FAIL) {
        fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
        return FAIL;
      } // end computefactor if
    } // end: if comoving coordinates

    GetUnits(&DensityUnits,&LengthUnits,&TemperatureUnits,&TimeUnits, &VelocityUnits,&MassUnits,Time);
    AccelUnits = LengthUnits/TimeUnits/TimeUnits;

    /* Loop over grid, adding acceleration to field. */
    for (dim = 0; dim < GridRank; dim++) {
      int n = 0;

      for (k = 0; k < GridDimension[2]; k++) {
        if (GridRank > 2)
          zpos=CellLeftEdge[2][k]+0.5*CellWidth[2][k]-DiskGravityPosition[2];
        if (dim == 2 && HydroMethod == Zeus_Hydro)
          zpos -= 0.5*CellWidth[2][k];

        for (j = 0; j < GridDimension[1]; j++) {
          if (GridRank > 1)
            ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j]-DiskGravityPosition[1];
          if (dim == 1 && HydroMethod == Zeus_Hydro)
            ypos -= 0.5*CellWidth[1][j];

          for (i = 0; i < GridDimension[0]; i++, n++) {
            xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]-DiskGravityPosition[0];
            if (dim == 0 && HydroMethod == Zeus_Hydro)
              xpos -= 0.5*CellWidth[0][i];

            /* Compute distance from center. */

            rsquared = xpos*xpos + ypos*ypos + zpos*zpos;

            double accelsph, accelcylR, accelcylz, zheight, xpos1, ypos1, zpos1;

            /* Compute z and r_perp (AngularMomentum is angular momentum 
             * and must have unit length). */

            /* magnitude of z = r.L in L direction */

            zheight=AngularMomentumx*xpos + AngularMomentumy*ypos + AngularMomentumz*zpos;

            /* position in plane of disk */

            xpos1=xpos-zheight*AngularMomentumx;
            ypos1=ypos-zheight*AngularMomentumy;
            zpos1=zpos-zheight*AngularMomentumz;

            radius = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1 + zheight*zheight);
            rcyl = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);
            radius = radius*LengthUnits;
            rcyl = rcyl*LengthUnits;
            accelsph = (GravConst)*MBulge*SolarMass/POW(radius+rBulge*Mpc,2)
                     + pi*GravConst*densDMConst*POW(rDMConst*Mpc,3)/POW(radius,2)
                       *(-2.0*atan(radius/rDMConst/Mpc)
                         +2.0*log(1.0+radius/rDMConst/Mpc)
                         +log(1.0+POW(radius/rDMConst/Mpc,2))
                        );
            accelcylR = GravConst*MSDisk*SolarMass*rcyl/sqrt(POW(POW(rcyl,2)
                      + POW(SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                      + POW(SDiskScaleHeightz*Mpc,2)),2),3));
            accelcylz = GravConst*MSDisk*SolarMass/sqrt(POW(zheight*LengthUnits,2)
                      + POW(SDiskScaleHeightz*Mpc,2))*zheight*LengthUnits/sqrt(POW(POW(rcyl,2)
                      + POW(SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                      + POW(SDiskScaleHeightz*Mpc,2)),2),3))
                        *(  SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                          + POW(SDiskScaleHeightz*Mpc,2))
                         )/AccelUnits;

             accelsph  = (radius ==0.0?0.0:fabs(accelsph )/(radius/LengthUnits)/AccelUnits);
             accelcylR = (rcyl   ==0.0?0.0:fabs(accelcylR)/(rcyl/LengthUnits)/AccelUnits);
             accelcylz = (zheight==0.0?0.0:fabs(accelcylz)*zheight/fabs(zheight));

             if (dim == 0)
               AccelerationField[0][n] -= (   accelsph*xpos
                                            + accelcylR*xpos1
                                            + accelcylz*AngularMomentumx);
             if (dim == 1)
               AccelerationField[1][n] -= (  accelsph*ypos
                                            + accelcylR*ypos1
                                            + accelcylz*AngularMomentumy);
             if (dim == 2)
               AccelerationField[2][n] -= (   accelsph*zpos
                                            + accelcylR*zpos1
                                            + accelcylz*AngularMomentumz);

          }
        }
      }  // end: loop over grid (i/j/k)
    } // end: loop over dims

    if (NumberOfParticles > 0 && GridRank != 3) {
        ENZO_FAIL("DiskGravity Requires 3D for use with particles");
    }

    if(NumberOfParticles > 0 && ParticleAcceleration[0] != NULL){
      for (int i = 0; i < NumberOfParticles; i++){

        // re-center coordinates relative to disk center, advancing vector by 1/2 timestep
        // these are actual (subgrid) positions
        // currently does not consider if mapping to grid cell pos is needed
        xpos = ParticlePosition[0][i] + 0.5*dtFixed*ParticleVelocity[0][i]/a - DiskGravityPosition[0];
        ypos = ParticlePosition[1][i] + 0.5*dtFixed*ParticleVelocity[1][i]/a - DiskGravityPosition[1];
        zpos = ParticlePosition[2][i] + 0.5*dtFixed*ParticleVelocity[2][i]/a - DiskGravityPosition[2];

        // model after grid loops
        rsquared = xpos*xpos + ypos*ypos + zpos*zpos;

        double accelsph, accelcylR, accelcylz, zheight, xpos1, ypos1, zpos1;

        // compute z and r_perp

        zheight = AngularMomentumx*xpos + AngularMomentumy*ypos + AngularMomentumz*zpos;

        // position in plane of disk
        xpos1 = xpos - zheight*AngularMomentumx;
        ypos1 = ypos - zheight*AngularMomentumy;
        zpos1 = zpos - zheight*AngularMomentumz;

        // again, copied from grid loops above
        radius = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1 + zheight*zheight);
        rcyl   = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);
        radius = radius*LengthUnits;
        rcyl   = rcyl*LengthUnits;

        accelsph = (GravConst)*MBulge*SolarMass/POW(radius+rBulge*Mpc,2)
                 + pi*GravConst*densDMConst*POW(rDMConst*Mpc,3)/POW(radius,2)
                 * (-2.0*atan(radius/rDMConst/Mpc)
                    +2.0*log(1.0+radius/rDMConst/Mpc)
                    +log(1.0+POW(radius/rDMConst/Mpc,2))
                   );
        accelcylR = GravConst*MSDisk*SolarMass*rcyl/sqrt(POW(POW(rcyl,2)
                  + POW(SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                  + POW(SDiskScaleHeightz*Mpc,2)),2),3));
        accelcylz = GravConst*MSDisk*SolarMass/sqrt(POW(zheight*LengthUnits,2)
                  + POW(SDiskScaleHeightz*Mpc,2))*zheight*LengthUnits/sqrt(POW(POW(rcyl,2)
                  + POW(SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                  + POW(SDiskScaleHeightz*Mpc,2)),2),3))
                    *(  SDiskScaleHeightR*Mpc+sqrt(POW(zheight*LengthUnits,2)
                      + POW(SDiskScaleHeightz*Mpc,2))
                     )/AccelUnits;

        accelsph  = (radius ==0.0?0.0:fabs(accelsph )/(radius/LengthUnits)/AccelUnits);
        accelcylR = (rcyl   ==0.0?0.0:fabs(accelcylR)/(rcyl/LengthUnits)/AccelUnits);
        accelcylz = (zheight==0.0?0.0:fabs(accelcylz)*zheight/fabs(zheight));


        ParticleAcceleration[0][i] -= (   accelsph*xpos
                                        + accelcylR*xpos1
                                        + accelcylz*AngularMomentumx);
        ParticleAcceleration[1][i] -= (   accelsph*ypos
                                        + accelcylR*ypos1
                                        + accelcylz*AngularMomentumy);
        ParticleAcceleration[2][i] -= (   accelsph*zpos
                                        + accelcylR*zpos1
                                        + accelcylz*AngularMomentumz);


      } // end: loop over particles
    } // end: check if particles + NULL

  } // end: if (DiskGravity)


  if (ExternalGravity) {

    if (PointSourceGravity)
      ENZO_FAIL("Cannot have both PointSourceGravity and ExternalGravity");
    /* Actually, you can, but it seems like a bad idea. 
       Correct if there is an exception to this */

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
    double Mvir = 4.0*M_PI*rhoc*POW(rs,3)*(log(1.0+c)-c/(1.0+c));
    
    float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, 
      TimeUnits = 1.0, VelocityUnits = 1.0;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	     &TimeUnits, &VelocityUnits, &MassUnits, Time);
    double AccelerationUnits = LengthUnits / POW(TimeUnits,2);
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
	    M = 4.0*M_PI*rhoc*POW(rs,3)*(log(1.0+x1)-x1/(1.0+x1));
	  }
	  else {
	    M = Mvir;
	  }
	  g = CGSGravConst*M/POW(r*LengthUnits,2);
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
	M = 4.0*M_PI*rhoc*POW(rs,3)*(log(1.0+x1)-x1/(1.0+x1));
      }
      else {
	M = Mvir;
      }
      g = CGSGravConst*M/POW(r*LengthUnits,2);
      g /= AccelerationUnits;

      ParticleAcceleration[0][i] += -g*xpos/r;
      ParticleAcceleration[1][i] += -g*ypos/r;
      ParticleAcceleration[2][i] += -g*zpos/r;
    }
    
  } // end if (ExternalGravity == 1)

  /* -----------------------------------------------------------------
     ExternalGravity > 9 : Acceleration from specified potential field.

     ExternalGravity == 10 : Flat rotation curve from Binney+Tremaine
     ExternalGravity == 20 : Point source gravity

     see Grid_AddExternalPotentialField.C for details.
     ----------------------------------------------------------------- */

  else if (ExternalGravity > 9) {

      /* Potential field size */
      int size1 = 1;
      for (dim = 0; dim < GridRank; dim++)
	size1 *= GravitatingMassFieldDimension[dim];

      /* Array for storying the external potential field */
      float *expotential = new float[size1];
      for (i = 0; i < size1; i++)
	expotential[i] = 0.0;

      /* Get external potential field for grid */
      if (this->AddExternalPotentialField(expotential) == FAIL) {
	fprintf(stderr, "Error in grid->AddExternalPotentialField().\n");
	return FAIL;
      }

      float *accel_field[MAX_DIMENSION];

      if (NumberOfBaryonFields > 0) {
    
	for (dim = 0; dim < GridRank; dim++)
	  accel_field[dim] = new float[size];   

	if (this->ComputeAccelerationsFromExternalPotential(
		  (HydroMethod == Zeus_Hydro) ? ZEUS_GRIDS : GRIDS, 
		  expotential, accel_field)  == FAIL) {
	  fprintf(stderr, "Error in grid->ComputeAccelerationForExternalPotential.\n");
	  return FAIL;
	}   

	for (dim = 0; dim < GridRank; dim++)
	  for (i = 0; i < size; i++)
	    AccelerationField[dim][i] += accel_field[dim][i];

	for (dim = 0; dim < GridRank; dim++) {
          delete [] accel_field[dim];
          accel_field[dim] = NULL;
	}

      } // end if (NumberOfBaryonFields)
 
      if (NumberOfParticles > 0) {

	for (dim = 0; dim < GridRank; dim++)
          accel_field[dim] = new float[NumberOfParticles];

	if (this->ComputeAccelerationsFromExternalPotential(PARTICLES, 
				  expotential, accel_field) == FAIL) {
	  fprintf(stderr, "Error in grid->ComputeAccelerationForExternalPotential.\n");
	  return FAIL;
	}  

	for (i = 0; i < NumberOfParticles; i++){
	  ParticleAcceleration[0][i] += accel_field[0][i];
	  ParticleAcceleration[1][i] += accel_field[1][i];
	  ParticleAcceleration[2][i] += accel_field[2][i];
	}

	for (dim = 0; dim < GridRank; dim++) {
          delete [] accel_field[dim];
          accel_field[dim] = NULL;
        }

      } // end if (NumberOfParticles)

      /* Cleanup */
      delete [] expotential;
      
   } // end if (ExternalGravity > 9)

  else {
      /* this is only reached if there are two types of external gravities - 
	 when you add a new one, this changes */
      ENZO_FAIL("should never get here! in Grid::ComputeAccelFieldExternal // ExternalGravity");
    }

  } // end if (ExternalGravity)


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
 
