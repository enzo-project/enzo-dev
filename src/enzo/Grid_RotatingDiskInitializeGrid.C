/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A ROTATING DISK SIMULATION)
/
/  written by: Elizabeth Tasker
/  date:       May 2012
/  modified1:  
/              
/  
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"
#include "units.h"

/* External routines */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

/* Internal routines */

float gas_velocity(FLOAT r, FLOAT cdensity, 
		  FLOAT totaldmmass, FLOAT rs, FLOAT hs, 
		  FLOAT dmc, FLOAT router);
	     
float average_density(FLOAT r, float cdensity, float rs, float hs, 
		      FLOAT cellwidth, FLOAT xpos, FLOAT ypos, FLOAT zpos);


int grid::RotatingDiskInitializeGrid(float RDScaleRadius,
				     float RDScaleHeight, 
				     float RDTemperature,
				     float RDDMConcentration, 
				     float RDTotalDMMass,
				     float RDCentralDensity,
				     float RDOuterRadius)
{

  
  /* create fields */
	
   int vel;
   NumberOfBaryonFields = 0;
   FieldType[NumberOfBaryonFields++] = Density;
   FieldType[NumberOfBaryonFields++] = TotalEnergy;
   if (DualEnergyFormalism)
     FieldType[NumberOfBaryonFields++] = InternalEnergy;
   vel = NumberOfBaryonFields;
   FieldType[NumberOfBaryonFields++] = Velocity1;
   if (GridRank > 1) 
     FieldType[NumberOfBaryonFields++] = Velocity2;
   if (GridRank > 2)
     FieldType[NumberOfBaryonFields++] = Velocity3;	


   /* Set various units. */

   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   GravitationalConstant = 4.0*pi*GravConst*MassUnits*pow(TimeUnits,2)/pow(LengthUnits,3);
  
   // NFW external potential required for this test problem
   PointSourceGravity = 2;

   /* Return if this grid is not being computed on this processor */
  
   if (MyProcessorNumber != ProcessorNumber) 
     return SUCCESS;
  
   /* declarations */
	
   int dim, i, j, k;
   FLOAT density, velocity[MAX_DIMENSION], temperature;
   FLOAT r, x, y = 0, z = 0, vx, vy, vz, velocity_mag;
   int n = 0;
   
		
   /* allocate memory for baryon fields */
		
   this->AllocateGrids();
   
   /* make disk position the same as external gravity center */

   for (dim = 0; dim < MAX_DIMENSION; dim++ ) 
     PointSourceGravityPosition[dim]= 0.5*(DomainRightEdge[dim]-DomainLeftEdge[dim]);

   
   /* Loop over the mesh. */
   
   for (k = 0; k < GridDimension[2]; k++)
     for (j = 0; j < GridDimension[1]; j++)
       for (i = 0; i < GridDimension[0]; i++, n++) {
	 
	 /* Compute physical position of each cell */
	
	 x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	 if (GridRank > 1)
	   y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	 if (GridRank > 2)
	   z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	 /* Background properties */
	 density = 1e-10;
	 temperature = RDTemperature;
	 for (dim = 0; dim < MAX_DIMENSION; dim++)
	   velocity[dim] = 0;
   
	 FLOAT xpos, ypos, zpos, drad;

	 /* Loop over dims if using Zeus (since vel's face-centered). */							
	 for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);  dim++) {
	 
	   /* Compute position from centre of galaxy */									
	   xpos = x-PointSourceGravityPosition[0] - 
	     (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	   ypos = y-PointSourceGravityPosition[1] -
	     (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	   zpos = z-PointSourceGravityPosition[2] -
	     (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

	   /* position in plane of disk */
	   drad = sqrt(xpos*xpos + ypos*ypos);

	   xpos = xpos/drad;
	   ypos = ypos/drad;
	   zpos = zpos/drad;

	   if (drad < RDOuterRadius ) 
	      if (dim == 0) { // calculate everything except velocities if using Zeus
		
		if (RDScaleHeight > CellWidth[0][0])
		  // calculate average density of cell in code units
		  density = average_density(drad, RDCentralDensity, RDScaleRadius, RDScaleHeight, 
				   CellWidth[0][0],  xpos*drad, ypos*drad, zpos*drad);
		else
		  density = RDCentralDensity*PEXP(-drad/RDScaleRadius)/POW(cosh(zpos*drad/CellWidth[0][0]), 2);

	      } // end if dim == 0
	  
	   velocity_mag = gas_velocity(drad, RDCentralDensity, RDTotalDMMass, RDScaleRadius, RDScaleHeight, RDDMConcentration, RDOuterRadius);
	    
	   /* Compute velocty: L x r_perp. */
						      
	   if (dim == 0 || dim == 1)
	     velocity[0] = -velocity_mag*ypos;
	   if (dim == 0 || dim == 2)
	     velocity[1] = velocity_mag*xpos;
	   if (dim == 0 || dim == 3)
	     velocity[2] = 0.0;
	   
	 } // end: loop over dims

	 BaryonField[0][n] = max(density, 1e-10);
								
	 /* Set energy (thermal and then total if necessary). */
	 BaryonField[1][n] = temperature/TemperatureUnits/((Gamma-1.0)*Mu);

	 if (DualEnergyFormalism)
	   BaryonField[2][n] = BaryonField[1][n];

	 if (HydroMethod != Zeus_Hydro)
	   for (dim = 0; dim < GridRank; dim++)
	     BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);
							
	 for (dim = 0; dim < GridRank; dim++)
	   BaryonField[vel+dim][n] = velocity[dim];

       } // end loop over grid
   
   return SUCCESS;
}


float average_density(FLOAT r, float cdensity, float rs, float hs, 
	     FLOAT cellwidth, FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
	// routine to return the average gas density in a grid cell
	// Routine samples density in r plane of grid and averages
	// Assumes all input units are code.
	
  int i,points;
  double den,r1,nx,ny;
	
  points = 100;
  den = cdensity*PEXP(-r/rs)/(POW(cosh(zpos/(hs)),2));
	
  for (i=0;i<points;i++)
    {
      nx = drand48()*cellwidth-cellwidth/2.0;
      ny = drand48()*cellwidth-cellwidth/2.0;
      r1 = sqrt(POW((xpos+nx),2)+POW((ypos+ny),2)); 
      den = den+cdensity*PEXP(-r1/rs)/(POW(cosh(zpos/(hs)),2));
    }
	
  double av_den = den/points;
	
  return av_den; //code units
	
}

float gas_velocity(FLOAT radius, FLOAT cdensity, FLOAT totaldmmass, FLOAT rscale, 
		   FLOAT hscale, FLOAT DMC, FLOAT router)
{

  /* calculates the circular velocity from mass of gas and dark matter */

  // cosmology values for dark matter halo
  double Redshift=1.0;	  
  double ExpansionFactor=1.0/(1.0+Redshift);
  double OmegaLambda=0.7;
  double OmegaMatter=0.3;
  double HubbleConstant=0.73;
  double OMEGA=OmegaLambda+OmegaMatter; // Flat Universe

  // units 
	
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1, time1=0.0;
  double MassUnits=1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, time1) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  double r_cgs = double(radius*LengthUnits);          // Radius [cm]
  double M_200 = double(totaldmmass*SolarMass);      // Virial Mass [g]
  
  // Hubble constant at z = Redshift
  double H = sqrt(HubbleConstant*100*HubbleConstant*100*(OmegaLambda+OmegaMatter*POW(ExpansionFactor,-3)-(OMEGA-1.)*POW(ExpansionFactor,-2)));     

  // Virial radius: equation from Navarro, Frank and White, 1997 gives r_200 in kpc. Convert to cm.
  double r_200 = (1.63e-2*POW(totaldmmass*MassUnits/SolarMass,1./3.)*POW((OmegaLambda+OmegaMatter*POW(ExpansionFactor, -3)-(OMEGA-1.0)*POW(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*POW(H,-2.0/3.0)*POW(100,2.0/3.0))*kpc;

  double massgas, massdm, masstotal, accel, vcirc;
  double fC = log(1.0+DMC)-DMC/(1.0+DMC);
  double rs = r_200/DMC;  //[m]

  if (radius < router) 
    massgas = 4.0*M_PI*double(hscale)*POW(double(rscale),2.0)*double(cdensity)*PEXP(-double(radius/rscale))*(PEXP(double(radius/rscale))-double(radius/rscale)-1.0)*MassUnits; 
  
  else // if we're past outer radius
    massgas = 4.0*M_PI*double(hscale)*POW(double(rscale),2.0)*double(cdensity)*PEXP(-double(router/rscale))*(PEXP(double(router/rscale))-double(router/rscale)-1.0)*MassUnits;

  massdm=(M_200/fC)*(log(1.0+r_cgs/rs)-(r_cgs/rs)/(1.0+r_cgs/rs));
  
  if (SelfGravity == 1)
    masstotal = massdm+massgas;
  else
    masstotal = massdm;

  PointSourceGravityConstant = FLOAT((M_200/fC)*(log(1.0+1.0)-1.0/(1.0+1.0)));
  PointSourceGravityCoreRadius = FLOAT(rs);

  accel = double(GravConst)*masstotal/(double(r_cgs)*double(r_cgs));
  vcirc = sqrt(r_cgs*accel);   
  
  return (FLOAT(vcirc/VelocityUnits));  //code units
   
}
