/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/  modified1:  Elizabeth Tasker, Oct, 2006 (tidied up)
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]


int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, float *MassUnits, FLOAT Time);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* Internal routines */

float gasdev();
float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, 
	     float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, 
	     float DMConcentration, FLOAT Time);
float av_den(FLOAT r, float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, 
	     FLOAT cellwidth, FLOAT z,FLOAT xpos, FLOAT ypos, FLOAT zpos);


static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits;

int grid::GalaxySimulationInitializeGrid(FLOAT DiskRadius,
					 float GalaxyMass,
					 float GasMass,
					 FLOAT DiskPosition[MAX_DIMENSION], 
					 FLOAT ScaleHeightz,
					 FLOAT ScaleHeightR, 
					 float DMConcentration,
					 float DiskTemperature,
					 float InitialTemperature,
					 float AngularMomentum[MAX_DIMENSION],
					 float UniformVelocity[MAX_DIMENSION], 
					 int UseMetallicityField, 
					 float GalaxySimulationInflowTime,
					 float GalaxySimulationInflowDensity,
					 int level)
{

 /* Return if this doesn't concern us. */

 if (ProcessorNumber != MyProcessorNumber) 
   return SUCCESS;

 /* declarations */

 int dim, i, j, k, m, field, disk, size, MetalNum, vel;
 float DiskDensity, DiskVelocityMag,dens2;

 // BWO: little indexing hack used near the end of this routine 
 if(DualEnergyFormalism){
   vel = 3;
 } else {
   vel = 2;
 }


 /* Set various units. */

 float CriticalDensity = 1, BoxLength = 1, mu = 0.6;
 FLOAT a, dadt, ExpansionFactor = 1;
 if (ComovingCoordinates) {
   CosmologyComputeExpansionFactor(Time, &a, &dadt);
   ExpansionFactor = a/(1.0+InitialRedshift);
   CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);
   CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
   BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
 }

 /* Set up inflow */
 if (GalaxySimulationInflowTime > 0.0){
   TimeActionType[0] = 2;
   TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
   TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
 }



 /* Set densities */

 float BaryonMeanDensity = 1.0;

 /* compute size of fields */

 size = 1;
 for (dim = 0; dim < GridRank; dim++)
   size *= GridDimension[dim];

 /* allocate fields */

 for (field = 0; field < NumberOfBaryonFields; field++)
   if (BaryonField[field] == NULL)
     BaryonField[field] = new float[size];

 /* Loop over the mesh. */

 float density, dens1, Velocity[MAX_DIMENSION],
   temperature, temp1, sigma, sigma1;
 FLOAT r, x, y = 0, z = 0;
 int n = 0;

 for (k = 0; k < GridDimension[2]; k++)
   for (j = 0; j < GridDimension[1]; j++)
     for (i = 0; i < GridDimension[0]; i++, n++) {

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	density = 1.0;
	temperature = temp1 = InitialTemperature;
	sigma = sigma1 = 0;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	/* Find distance from center. */

	r = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		 POW(fabs(y-DiskPosition[1]), 2) +
		 POW(fabs(z-DiskPosition[2]), 2) );
	r = max(r, 0.1*CellWidth[0][0]);

	if (r < DiskRadius) {

	  FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;

	  /* Loop over dims if using Zeus (since vel's face-centered). */

	  for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
	       dim++) {

	    /* Compute position. */

	    xpos = x-DiskPosition[0] - 
	      (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	    ypos = y-DiskPosition[1] -
	      (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	    zpos = z-DiskPosition[2] -
	      (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);
	    
	    /* Compute z and r_perp (AngularMomentum is angular momentum 
	       and must have unit length). */    

	    /* magnitude of z = r.L in L direction */

	    zheight = AngularMomentum[0]*xpos + 
	              AngularMomentum[1]*ypos +
	              AngularMomentum[2]*zpos;

	    /* position in plane of disk */

	    xpos1 = xpos - zheight*AngularMomentum[0];
	    ypos1 = ypos - zheight*AngularMomentum[1];
	    zpos1 = zpos - zheight*AngularMomentum[2];
	    drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);


	    /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	    xpos1 = xpos1/drad;
	    ypos1 = ypos1/drad;
	    zpos1 = zpos1/drad;

	    /* If we're above the disk, then exit. */
	    
	    DiskDensity = (GasMass*SolarMass/(8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))/DensityUnits;   //Code units (rho_0)

	    DiskVelocityMag = gasvel(drad, DiskDensity, ExpansionFactor, GalaxyMass, ScaleHeightR, ScaleHeightz, DMConcentration, Time);


	    if (dim == 0)

	      if (ScaleHeightz*Mpc/LengthUnits > CellWidth[0][0])
		{

		  dens1=av_den(drad*LengthUnits, DiskDensity*DensityUnits, ScaleHeightR*Mpc, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits, zheight*LengthUnits,xpos1*drad*LengthUnits,ypos1*drad*LengthUnits,zpos1*drad*LengthUnits);

		}
	      else
		{
		  dens1 = DiskDensity*exp(-drad/(ScaleHeightR*Mpc/LengthUnits))/POW(cosh(zheight/CellWidth[0][0]), 2);
		}	  
	    dens2 = DiskDensity*exp(-drad/(ScaleHeightR*Mpc/LengthUnits))/POW(cosh(zheight/max((ScaleHeightz*Mpc/LengthUnits), CellWidth[0][0])), 2);
	    
	    if (dens2 < density)
	      break;

	    /* Compute velocity magnitude (divided by drad). 
	       This assumes PointSourceGravityPosition and Disk center 
	       are the same. */

	    /* Compute velocty: L x r_perp. */

	    if (dim == 0 || dim == 1)
	      Velocity[0] = DiskVelocityMag*(AngularMomentum[1]*zpos1 -
					     AngularMomentum[2]*ypos1);
	    if (dim == 0 || dim == 2)
	      Velocity[1] = DiskVelocityMag*(AngularMomentum[2]*xpos1 -
					     AngularMomentum[0]*zpos1);
	    if (dim == 0 || dim == 3)
	      Velocity[2] = DiskVelocityMag*(AngularMomentum[0]*ypos1 -
					     AngularMomentum[1]*xpos1);
	    
	  } // end: loop over dims

	   	    
	    /* If the density is larger than the background (or the previous
	       disk), then set the velocity. */

	  if (dens2 > density) {
	    density = dens1;
	    if (temp1 == InitialTemperature)
	      temp1 = DiskTemperature;
	    temperature = temp1;
	    sigma = sigma1;
	  }

	} // end: if (r < DiskRadius)
	
	/* Set density. */

	BaryonField[0][n] = density*BaryonMeanDensity;
	
	if (UseMetallicityField)
	  for (i = 0; i < size; i++)
	    BaryonField[MetalNum][i] = 1.0e-10;

	/* Set Velocities. */

	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits/
                           ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];
	
	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);

	if (BaryonField[1][n] <= 0)
	  printf("n = %d  temp = %g   e = %g\n", 0, temperature, 
	       BaryonField[1][0]);


     } // end loop over grid

 return SUCCESS;

}


float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, float DMConcentration, FLOAT Time)
{

 double OMEGA=OmegaLambdaNow+OmegaMatterNow;                 //Flat Universe

 double r = radius*LengthUnits/100;    // Radius [m]

 double M_200 = GalaxyMass*SolarMass/1000.0;      // Virial Mass [kg]

 double H = sqrt(HubbleConstantNow*100*HubbleConstantNow*100*(OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor,-3)-(OMEGA-1.)*POW(ExpansionFactor,-2)));                                

 double r_200 = (1.63e-2*POW(GalaxyMass,1.0/3.0)*POW((OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor, -3)-(OMEGA-1.0)*POW(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*POW(H,-2.0/3.0)*POW(100,2.0/3.0))*Mpc/1.0e5;
 //virial radius [m]: M_200/M_Solar = GalaxyMass

 double M_gas, M_DM, M_Tot, Acc, V_Circ;
 double f_C = log(1.0+DMConcentration)-DMConcentration/(1.0+DMConcentration);
 double r_s = r_200/DMConcentration;  //[m]

 // Mass of gas disk and DM at given radius

     M_gas=8.0*M_PI*ScaleHeightz*Mpc/100*ScaleHeightR*Mpc/100*ScaleHeightR*Mpc/100*DiskDensity*DensityUnits*1000*exp(-r/(ScaleHeightR*Mpc/100))*(exp(r/(ScaleHeightR*Mpc/100))-r/(ScaleHeightR*Mpc/100)-1.0);

     M_DM=(M_200/f_C)*(log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s));

     if (SelfGravity==1){
	M_Tot=M_DM+M_gas;
     }
     else{
	M_Tot=M_DM;
     }

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1, MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  double MassUnitsDouble=1.0;

  if(ComovingCoordinates)
    MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);

  // Set the point source gravity parameters.  This is the DM mass (in g)
  //   within rs.  The core radius to rs in cm.
  //
  // BWO 10 July 2009: Both of these values are now converted to code units, because 
  // otherwise the values go over 32-bit precision.  This is used in
  // Grid::ComputeAccelerationFieldExternal, and converted back to CGS where needed.
  //

  PointSourceGravityConstant = (M_200/f_C)*(log(1.0+1.0)-1.0/(1.0+1.0))*1000.0 / MassUnitsDouble;
  PointSourceGravityCoreRadius = r_s*100.0 / LengthUnits;

  /*
  fprintf(stderr,"Grid::GalaxySimulationInitializeGrid:  %d  %e  %e\n",MyProcessorNumber,MassUnitsDouble, LengthUnits);
  fprintf(stderr,"  PointSourceGravityConstant = %e  %d\n",PointSourceGravityConstant,MyProcessorNumber);
  fprintf(stderr,"  PointSourceGravityCoreRadius = %e  %d\n",PointSourceGravityCoreRadius,MyProcessorNumber);
  */

 // Force per unit mass on disk (i.e. acceleration) [ms-2]

     Acc=((GravConst/1000.0)*M_Tot)/(r*r);

 // Magnitude of Circular Velocity of disk 

     V_Circ = sqrt(r*Acc)*100;       //cms-1

     /*      printf("r = %g  M_Tot = %g  Acc = %g  M_DM = %g  M_gas = %g  f_C = %g\n",
	     r, M_Tot, Acc, M_DM, M_gas, f_C);
     printf("r_s = %g  DMConcentration = %g  r_200 = %g  r/r_s = %g\n",
	     r_s, DMConcentration, r_200, r/r_s);
     printf("EF = %g  H = %g  OMEGA = %g\n", ExpansionFactor, H, OMEGA);
     printf("radius = %g  v_circ = %g\n", radius, V_Circ);  */

     return (V_Circ/VelocityUnits);  //code units
}

float av_den(FLOAT r, float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
 // routine to return the average gas density in a grid cell
 // Routine samples density in r plane of grid and averages
 // Assumes all input units are CGS!!

 int i,points;
 double den,r1,nx,ny,nz;

 points = 100;
 den = DiskDensity*exp(-r/ScaleHeightR)/(POW(cosh(z/(2.0*ScaleHeightz)),2));

 for (i=0;i<points;i++)
   {
     nx = drand48()*cellwidth-cellwidth/2.0;
     ny = drand48()*cellwidth-cellwidth/2.0;
     nz = drand48()*cellwidth-cellwidth/2.0;
     r1 = sqrt(POW((xpos+nx),2)+POW((ypos+ny),2)+POW((zpos+nz),2)); 
     den = den+DiskDensity*exp(-r1/ScaleHeightR)/(POW(cosh(z/(2.0*ScaleHeightz)),2));
   }

 double av_den = den/points;

 return (av_den/DensityUnits); //code unites

}
