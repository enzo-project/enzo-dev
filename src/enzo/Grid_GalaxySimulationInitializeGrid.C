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
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
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
	     float *VelocityUnits, FLOAT Time);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* Internal routines */

float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, 
	     float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, 
	     float DMConcentration, FLOAT Time);
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth);
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);

/* Internal Routines for Disk Potential Setup */
float HaloGasDensity(FLOAT);
float HaloGasTemperature(FLOAT);
double DiskPotentialCircularVelocity(FLOAT cellwidth,FLOAT z,FLOAT density,FLOAT &temperature);
double trapzd(double (func)(), double a, double b, int n);
double qromb(double (*func)(double), double a, double b);
void polint(double xa[],double ya[],int n,double x,double *y,double *dy);
double func1(double zint);
double func2(double zint);
double func3(double zint);
double func4(double zint);
double *vector(int nl,int nh);
void free_vector(double *v,int nl,int nh);
static double drcyl;
static double r2;

static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits, MassUnits;

double gScaleHeightR, gScaleHeightz, densicm, MgasScale, Picm, TruncRadius, SmoothRadius, SmoothLength,Ticm;
double GalaxySimulationGasHalo, GalaxySimulationGasHaloScaleRadius, GalaxySimulationGasHaloDensity;

int grid::GalaxySimulationInitializeGrid(FLOAT DiskRadius,
					 float GalaxyMass,
					 float GasMass,
					 FLOAT DiskPosition[MAX_DIMENSION], 
					 FLOAT ScaleHeightz,
					 FLOAT ScaleHeightR,
					 FLOAT GalaxyTruncationRadius, 
					 float DMConcentration,
					 float DiskTemperature,
					 float InitialTemperature,
					 float UniformDensity,
					 int   GasHalo,
					 float GasHaloScaleRadius,
					 float GasHaloDensity,
					 float AngularMomentum[MAX_DIMENSION],
					 float UniformVelocity[MAX_DIMENSION], 
					 int UseMetallicityField, 
					 float GalaxySimulationInflowTime,
					 float GalaxySimulationInflowDensity,
					 int level,
					 float GalaxySimulationCR )
{
 /* declarations */

  int dim, i, j, k, m, field, disk, size, MetalNum, MetalIaNum, vel;
 int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
   DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
 float DiskDensity, DiskVelocityMag;
  int CRNum, DensNum;

  /* global-scope variables for disk potential functions (would be better if not global) */

  gScaleHeightR = ScaleHeightR;
  gScaleHeightz = ScaleHeightz;
  densicm = UniformDensity;
  MgasScale = GasMass;
	Ticm = InitialTemperature;
  Picm = kboltz*UniformDensity*Ticm/(0.6*mh);
  TruncRadius = GalaxyTruncationRadius;
  SmoothRadius = TruncRadius*.02/.026;
  SmoothLength = TruncRadius - SmoothRadius;

	GalaxySimulationGasHalo = GasHalo;
	GalaxySimulationGasHaloScaleRadius = GasHaloScaleRadius;
	GalaxySimulationGasHaloDensity = GasHaloDensity;

  /* create fields */

  NumberOfBaryonFields = 0;
  DensNum = NumberOfBaryonFields;
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
  if (HydroMethod == MHD_RK) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
    FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    if (UseDivergenceCleaning) {
      FieldType[NumberOfBaryonFields++] = Phi_pField;
    }
  }

  /* If cosmic rays present, set up field */
  CRNum = NumberOfBaryonFields;
  if( CRModel )
    FieldType[NumberOfBaryonFields++] = CRDensity;

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  if (UseMetallicityField)
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* fake it with metals */
  if (StarMakerTypeIaSNe)
    FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

 /* Return if this doesn't concern us. */

 if (ProcessorNumber != MyProcessorNumber) 
   return SUCCESS;

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
 } else if( PointSourceGravity ){
   ENZO_FAIL("ERROR IN GALAXY SIM GRID INITIALIZE: non-cosmology units not supported for point source gravity");
 } else {
   if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
  } // end get units error if  
 } // end units if/else

	/* correct background density if it's not given in code units */
	if( UniformDensity < 1.0E-10 ){
		UniformDensity /= DensityUnits;
		if( debug && MyProcessorNumber == ROOT_PROCESSOR ) 
			fprintf(stdout,"Converting GalaxySimulationUniformDensity = %"GSYM" from CGS to code units\n",UniformDensity);
	} // end uniform density if

 /* Set up inflow */
 if (GalaxySimulationInflowTime > 0.0){
   TimeActionType[0] = 2;
   TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
   TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
 }

 /* compute size of fields */

 size = 1;
 for (dim = 0; dim < GridRank; dim++)
   size *= GridDimension[dim];

 /* allocate fields */

 for (field = 0; field < NumberOfBaryonFields; field++)
   if (BaryonField[field] == NULL)
     BaryonField[field] = new float[size];

 /* set metals to small value */

  if (UseMetallicityField)
    for (i = 0; i < size; i++)
      BaryonField[MetalNum][i] = 1.0e-10;

 /* Loop over the mesh. */

 float density, dens1, Velocity[MAX_DIMENSION];
 FLOAT temperature, temp1, init_temp;
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

	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	/* Find distance from center. */

	r = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		 POW(fabs(y-DiskPosition[1]), 2) +
		 POW(fabs(z-DiskPosition[2]), 2) );
	r = max(r, 0.1*CellWidth[0][0]);

	density = HaloGasDensity(r)/DensityUnits;
	temperature = temp1 = init_temp = HaloGasTemperature(r);

	if (r < DiskRadius) {

	  FLOAT xpos, ypos, zpos, zheight, drad; 
	  float CellMass;
	  FLOAT xhat[3];
	  FLOAT yhat[3];

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

	    xhat[0] = xpos - zheight*AngularMomentum[0];
	    xhat[1] = ypos - zheight*AngularMomentum[1];
	    xhat[2] = zpos - zheight*AngularMomentum[2];
	    drad = sqrt(xhat[0]*xhat[0] + xhat[1]*xhat[1] + xhat[2]*xhat[2]);
      drcyl = drad;

	    /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	    xhat[0] = xhat[0]/drad;
	    xhat[1] = xhat[1]/drad;
	    xhat[2] = xhat[2]/drad;

	    /* Find another vector perpendicular to r_perp and AngularMomentum */

	    yhat[0] = AngularMomentum[1]*xhat[2] - AngularMomentum[2]*xhat[1];
	    yhat[1] = AngularMomentum[2]*xhat[0] - AngularMomentum[0]*xhat[2];
	    yhat[2] = AngularMomentum[0]*xhat[1] - AngularMomentum[1]*xhat[0];

	    /* generate rotation matrix */
	    FLOAT inv[3][3],temp;
	    int i,j;
	    
	    // matrix of basis vectors in coordinate system defined by the galaxy
	    inv[0][0] = xhat[0]; inv[0][1] = yhat[0]; inv[0][2] = AngularMomentum[0];
	    inv[1][0] = xhat[1]; inv[1][1] = yhat[1]; inv[1][2] = AngularMomentum[1];
	    inv[2][0] = xhat[2]; inv[2][1] = yhat[2]; inv[2][2] = AngularMomentum[2];
	    
	    // Matrix is orthogonal by construction so inverse = transpose
	    for (i=0;i<3;i++)
	      for (j=i+1;j<3;j++)
		{
		  temp = inv[i][j];
		  inv[i][j] = inv[j][i];
		  inv[j][i] = temp;
		}

	    if( fabs(drcyl*LengthUnits/Mpc) > TruncRadius ){
	      dens1 = 0.0;
	      break;
	    }
		  
	    DiskDensity = (GasMass*SolarMass/(8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))/DensityUnits;   //Code units (rho_0) DENSITY NEEDS CUTOFF

	    if (PointSourceGravity > 0 )
	      DiskVelocityMag = gasvel(drad, DiskDensity, ExpansionFactor, GalaxyMass, ScaleHeightR, ScaleHeightz, DMConcentration, Time);
	    else if( DiskGravity > 0 ){
	      CellMass = gauss_mass(drad*LengthUnits,zheight*LengthUnits, xpos*LengthUnits, ypos*LengthUnits, zpos*LengthUnits, inv,
				    DiskDensity*DensityUnits,ScaleHeightR*Mpc, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits);

	      dens1 = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;

	      DiskVelocityMag = DiskPotentialCircularVelocity(CellWidth[0][0], zheight*LengthUnits, dens1, temp1);
	    }
	    if (PointSourceGravity*DiskGravity != FALSE ) 
	      ENZO_FAIL("Cannot activate both PointSource and Disk gravity options for Isolated Galaxy");

	    if (dim == 0) {
	      CellMass = gauss_mass(drad*LengthUnits,zheight*LengthUnits, xpos*LengthUnits, ypos*LengthUnits, zpos*LengthUnits, inv, 
				    DiskDensity*DensityUnits,ScaleHeightR*Mpc, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits);
	      dens1 = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;
	    }

	    /* If we're above the disk, then exit. */

	    if (dens1 < density)
	      break;

	    /* Compute velocity magnitude (divided by drad). 
	       This assumes PointSourceGravityPosition and Disk center 
	       are the same. */

	    /* Compute velocty: L x r_perp. */

	    if (dim == 0 || dim == 1)
	      Velocity[0] = DiskVelocityMag*(AngularMomentum[1]*xhat[2] -
					     AngularMomentum[2]*xhat[1]);
	    if (dim == 0 || dim == 2)
	      Velocity[1] = DiskVelocityMag*(AngularMomentum[2]*xhat[0] -
					     AngularMomentum[0]*xhat[2]);
	    if (dim == 0 || dim == 3)
	      Velocity[2] = DiskVelocityMag*(AngularMomentum[0]*xhat[1] -
					     AngularMomentum[1]*xhat[0]);
	    
	  } // end: loop over dims

	   	    
	    /* If the density is larger than the background (or the previous
	       disk), then set the velocity. */

	  if (dens1 > density && fabs(drcyl*LengthUnits/Mpc) <= TruncRadius ) {
	    density = dens1;
	    if (temp1 == init_temp)
	      temp1 = DiskTemperature;
	    temperature = temp1;
	    if( temperature > 1e7 )
	      temperature = init_temp;
	    if( UseMetallicityField ) // This should be converted to a general color field at some point - this obviously breaks metallicity feature
	      BaryonField[MetalNum][n] = density;
	  }

	} // end: if (r < DiskRadius)
	
	/* Set density. */

	BaryonField[0][n] = density;
	
	if (StarMakerTypeIaSNe)
	  for (i = 0; i < size; i++)
	    BaryonField[MetalIaNum][i] = 1.0e-10;

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

     if( CRModel )
       BaryonField[CRNum][n] = BaryonField[DensNum][n] * GalaxySimulationCR;



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

     M_gas=8.0*M_PI*ScaleHeightz*Mpc/100*ScaleHeightR*Mpc/100*ScaleHeightR*Mpc/100*DiskDensity*DensityUnits*1000*PEXP(-r/(ScaleHeightR*Mpc/100))*(PEXP(r/(ScaleHeightR*Mpc/100))-r/(ScaleHeightR*Mpc/100)-1.0);

     M_DM=(M_200/f_C)*(log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s));

     if (SelfGravity==1){
	M_Tot=M_DM+M_gas;
     }
     else{
	M_Tot=M_DM;
     }

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
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


// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth)
{
  
  FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
  FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
  FLOAT xResult [5];
  FLOAT yResult [5];
  float Mass = 0;
  FLOAT xrot,yrot,zrot;
  int i,j,k;
  FLOAT rrot;
	
  for (i=0;i<5;i++) {

      xResult[i] = 0.0;
      for (j=0;j<5;j++) {

	  yResult[j] = 0.0;
	  for (k=0;k<5;k++) {

	      rot_to_disk(xpos+EvaluationPoints[i]*cellwidth/2.0,ypos+EvaluationPoints[j]*cellwidth/2.0,zpos+EvaluationPoints[k]*cellwidth/2.0,xrot,yrot,zrot,inv);
	      rrot = sqrt(POW(xrot,2)+POW(yrot,2));

	      if( PointSourceGravity > 0 )
		yResult[j] += cellwidth/2.0*Weights[k]*PEXP(-rrot/ScaleHeightR)/POW(cosh(zrot/(2.0*ScaleHeightz)),2);
	      else if( DiskGravity > 0 ){
		if( rrot/Mpc < SmoothRadius )
		  yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz);
		else if( rrot/Mpc < TruncRadius )
		  yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz)*0.5*(1.0+cos(pi*(rrot-SmoothRadius*Mpc)/(SmoothLength*Mpc)));
	      } // end disk gravity if

	  }
	  xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
      }
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
  }  
  Mass *= DiskDensity;
  return Mass;
}

//Finds coordinates in rotated coordinate system
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3])
{
  xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
  yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
  zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];
}


double DiskPotentialDarkMatterMass(FLOAT R){
/*
 * 	computes dark matter mass enclosed within spherical radius R
 *	for potential in Mori & Burkert 2000, consistent with eq
 *
 *		rho = rho0 * r0**3 / ( (r + r0)*(r**2 + r0**2 ) )
 *			
 *	Parameters:
 *	-----------
 *		R - Spherical radius (code units)
 *
 * 	Returns: Mass, in grams
 */
	FLOAT R0 = DiskGravityDarkMatterR*Mpc,x=R/R0*LengthUnits;
	double M0 = pi*DiskGravityDarkMatterDensity*R0*R0*R0;

	return M0*(-2.0*atan(x)+2.0*log(1+x)+log(1.0+x*x));
} // end DiskPotentialDarkMatterMass


float HaloGasTemperature(FLOAT R){
/*
 *	computes halo temperature, assuming gas particles follow
 *	KE = 1/2 PE assuming DM potential given in 
 *	DiskPotentialDarkMatterMass() above. 
 *
 *	Parameters:
 *	-----------
 *		R - spherical radius (code units)
 *
 *	Returns: Temperature, Kelvin
 */
	if(GalaxySimulationGasHalo)
		return GravConst*DiskPotentialDarkMatterMass(R)*0.6*mh/(3.0*kboltz*R*LengthUnits);
	return Ticm;
}


float DiskPotentialGasDensity(FLOAT r,FLOAT z){
/*
 *	computes gas density within galaxy disk, according to eq
 *
 * 		(Mgas/8*pi*a^2*b)*sech(r/a)*sech*(z/b)
 *
 * 	Smoothed by a cosine fcn beyond SmoothRadius
 *
 * 	Parameteres:
 * 	------------
 * 		r - cylindrical radius (code units)
 * 		z - cylindrical height (code units)
 *
 * 	Returns: density (in grams/cm^3)
 *
 */
	double density = MgasScale*SolarMass/(8.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc);
	density /= (cosh(r*LengthUnits/gScaleHeightR/Mpc)*cosh(z*LengthUnits/gScaleHeightz/Mpc));

	if(fabs(r*LengthUnits/Mpc) > SmoothRadius && fabs(r*LengthUnits/Mpc) <= TruncRadius)
		density *= 0.5*(1.0+cos(pi*(r*LengthUnits-SmoothRadius*Mpc)/(SmoothLength*Mpc)));
	return density;
} // end DiskPotentialGasDensity


float HaloGasDensity(FLOAT R){
/*
 *	computes density of gaseous halo, assuming hydro equilibrium
 *	and temperature profile given above in HaloGasTemperature().
 *
 * 	Paramaters:
 * 	-----------
 * 		R - spherical radius, code units
 *
 * 	Returns: density, grams/cm^3
 */
	if(GalaxySimulationGasHalo){
		double T0,haloDensity;
		T0 = HaloGasTemperature(GalaxySimulationGasHaloScaleRadius*Mpc/LengthUnits);
		haloDensity = GalaxySimulationGasHaloDensity*(T0/HaloGasTemperature(R));
		haloDensity /= pow((R*LengthUnits/GalaxySimulationGasHaloScaleRadius/Mpc),3);
		return min(haloDensity,GalaxySimulationGasHaloDensity);
	}
	return densicm;
} // end HaloGasDensity


double findZicm(FLOAT r){
  /*  
   *  Finds the height above the disk plane where the disk gas density
   *  matches the halo's gas density (using bisection)
   *
   *  Parameters:
   *  -----------
   *  	r - cylindrical radius (code units)
   *
   *  Returns: zicm, edge of disk, (code units)
   */

  static const double X_TOL = 1e-7*Mpc/LengthUnits; // sub pc resolution
  static const int MAX_ITERS = 50; int iters=0;

  double z_lo = 0.0,z_hi = 0.01*Mpc/LengthUnits,z_new,f_lo,f_hi,f_new;
  f_hi = DiskPotentialGasDensity(r,z_hi) - HaloGasDensity(sqrt(r*r+z_hi*z_hi)); // -ve
  f_lo = DiskPotentialGasDensity(r,z_lo) - HaloGasDensity(sqrt(r*r+z_lo*z_lo)); // +ve

  if(f_lo < 0.0) return 0.0; // beyond the disk
  if(f_hi > 0.0) ENZO_FAIL("ERROR IN GALAXY INITIALIZE: HALO IS UNDER-PRESSURIZED");

  while(iters++ < MAX_ITERS ){

    z_new = (z_hi+z_lo)/2.0;
    f_new = DiskPotentialGasDensity(r,z_new)
            - HaloGasDensity(sqrt(r*r+z_new*z_new));

    if( fabs(f_new) == 0.0 ) return z_new;
    if( f_new*f_lo > 0.0 ){
      z_lo = z_new; f_lo = f_new;
    }
    else{
      z_hi = z_new; f_hi = f_new;
    }
    if( fabs(z_hi - z_lo) <= X_TOL ) return z_new;
  }

  ENZO_FAIL("ERROR IN GALAXY INITIALIZE: findZicm FAILED TO CONVERGE");
  return -1.0;
}


/* 
 *	DISK POTENTIAL CIRCULAR VELOCITY
 */
float DiskPotentialCircularVelocity(FLOAT cellwidth, FLOAT z, FLOAT density, FLOAT &temperature)
{

	extern double drcyl;
	double func1(double zint);       //(density times Stellar bulge force)
	double func2(double zint);     //(density times stellar disk force)
	double func3(double zint);       //func1 but for r2
	double func4(double zint);      //func2 but for r2

	double Pressure,Pressure2,zicm,zicm2,zicmf=0.0,zsmall=0.0,
		zicm2f=0.0,zint,FdPdR,FtotR,denuse,rsph,vrot,bulgeComp,rsph_icm;

	r2=(drcyl+0.01*cellwidth)*LengthUnits;
	rsph=sqrt(pow(drcyl*LengthUnits,2)+pow(z,2));

	/*	Determine zicm: the height above the disk where rho -> rho_ICM,
	 *	use this to find P_icm and dP_icm  */
	if (fabs(drcyl*LengthUnits/Mpc) <= SmoothRadius) {

		zicm  = findZicm(drcyl)*LengthUnits;
		zicm2 = findZicm(r2/LengthUnits)*LengthUnits;

		if( fabs(z) < fabs(zicm) ){
			bulgeComp = (DiskGravityStellarBulgeMass==0.0?0.0:qromb(func1, fabs(zicm), fabs(z)));
			Pressure= bulgeComp + qromb(func2, fabs(zicm), fabs(z));
			bulgeComp = (DiskGravityStellarBulgeMass==0.0?0.0:qromb(func3, fabs(zicm2), fabs(z)));
			Pressure2= bulgeComp + qromb(func4, fabs(zicm2), fabs(z));
		}  // end |z| < |zicm| if
	}  else {
    if (fabs(drcyl*LengthUnits/Mpc) <= TruncRadius ) {

			zicm  = findZicm(drcyl)*LengthUnits;
			zicm2 = findZicm(r2/LengthUnits)*LengthUnits;

/*
			if ( HaloGasDensity(sqrt(drcyl*drcyl+z*z)) >= DiskPotentialGasDensity(drcyl,z)
					&& fabs(z) < zicm) {
				printf("small density zicm = %g, z = %g\n", zicm/Mpc, z/Mpc);
			} // end small density if
*/ // FIXME

			if (fabs(z) < fabs(zicm)) {

				bulgeComp = (DiskGravityStellarBulgeMass==0.0?0.0:qromb(func1, fabs(zicm), fabs(z)));
				Pressure  = (bulgeComp+ qromb(func2, fabs(zicm), fabs(z)))
				            *(0.5*(1.0+cos(pi*(drcyl*LengthUnits-SmoothRadius*Mpc)/(SmoothLength*Mpc))));
				bulgeComp = (DiskGravityStellarBulgeMass==0.0?0.0:qromb(func3, fabs(zicm2), fabs(z)));
    		Pressure2 = (bulgeComp + qromb(func4, fabs(zicm2), fabs(z)))
				            *(0.5*(1.0+cos(pi*(r2-SmoothRadius*Mpc)/(SmoothLength*Mpc))));
			} // end |z| < |zicm| if

  	} // end r_cyle < TruncRadius if
	} // end r_cyl < SmoothRadius if/else

	denuse = density*DensityUnits; 
	if (Pressure < 0.0 && fabs(drcyl)*LengthUnits/Mpc <= TruncRadius && fabs(z) <= fabs(zicm)) {
		fprintf(stderr,"neg pressure:  P = %"FSYM", z = %"FSYM", r = %"FSYM"\n", Pressure, z/Mpc, drcyl*LengthUnits/Mpc);
	}
	if (fabs(drcyl)*LengthUnits/Mpc >= TruncRadius || fabs(zicm) <= fabs(z)){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = HaloGasDensity(rsph);
	}
	if (Pressure2 <= 0.0 && Pressure <= 0.0){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = HaloGasDensity(rsph);
	}
	if (Pressure <= 0.0) {
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = HaloGasDensity(rsph);
	}
	if (denuse < HaloGasDensity(rsph)) {
		fprintf(stderr,"denuse small:  %"FSYM"\n", denuse);
	}
	rsph_icm = sqrt(drcyl*drcyl+pow(zicm/LengthUnits,2));
	Picm = HaloGasDensity(rsph_icm)*kboltz*HaloGasTemperature(rsph_icm)/(0.6*mh);
	temperature=0.6*mh*(Picm+Pressure)/(kboltz*denuse);

	/* Calculate pressure gradient */
	FdPdR = (Pressure2 - Pressure)/(r2-drcyl*LengthUnits)/density; 

	/* Calculate Gravity = Fg_DM + Fg_StellarDisk + Fg_StellaDiskGravityStellarBulgeR */
	FtotR  = (-pi)*GravConst*DiskGravityDarkMatterDensity*pow(DiskGravityDarkMatterR*Mpc,3)/pow(rsph,3)*drcyl*LengthUnits
	         *(-2.0*atan(rsph/DiskGravityDarkMatterR/Mpc)+2.0*log(1.0+rsph/DiskGravityDarkMatterR/Mpc)
	           +log(1.0+pow(rsph/DiskGravityDarkMatterR/Mpc,2)));
	FtotR += -GravConst*DiskGravityStellarDiskMass*SolarMass*drcyl*LengthUnits
             /sqrt(pow(pow(drcyl*LengthUnits,2)+pow(DiskGravityStellarDiskScaleHeightR*Mpc
	                   +sqrt(pow(z,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)),2),3));
	FtotR += -GravConst*DiskGravityStellarBulgeMass*SolarMass
	           /pow(sqrt(pow(z,2)+pow(drcyl*LengthUnits,2))+DiskGravityStellarBulgeR*Mpc,2)*drcyl*LengthUnits/sqrt(pow(z,2)
	                +pow(drcyl*LengthUnits,2));


	if (temperature < 0.0) fprintf(stderr,"temp = %g, P = %g, z = %g, zicm = %g, zicmf=%g, zsmall=%g, drcyl = %g\n", 
		temperature, Pressure, z/Mpc, zicm/Mpc, zicmf, zsmall, drcyl*LengthUnits/Mpc);
	if ((FtotR - FdPdR) > 0.0) { 
		fprintf(stderr,"FtotR = %g, FdPdR = %g, P = %g,P2 = %g, Picm = %g, dr = %g, drcyl = %g, z = %g\n", 
			FtotR, FdPdR, Pressure, Pressure2, Picm, r2-drcyl*LengthUnits, drcyl*LengthUnits/Mpc, z/Mpc);
   	FdPdR = 0.0;
	} // end FtotR - FdPdr > 0.0 if

	/* Find circular velocity by balancing FG and dPdR against centrifugal force */
	vrot=sqrt(-drcyl*LengthUnits*(FtotR-FdPdR));
	if ((denuse == densicm)) vrot = 0.0;

  return (vrot/VelocityUnits); //code units
} // end DiskPotentialCircularVelocity


// the two functions integrated by qromb

double func1(double zint)
{

  extern double drcyl;
  return (-MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(zint)/gScaleHeightz/Mpc)*GravConst*DiskGravityStellarBulgeMass*SolarMass/pow((sqrt(pow(zint,2)+pow(drcyl*LengthUnits,2))+DiskGravityStellarBulgeR*Mpc),2)*fabs(zint)/sqrt(pow(zint,2)+pow(drcyl*LengthUnits,2)));
}

double func2(double zint)
{
  extern double drcyl;
  return (-MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(zint)/gScaleHeightz/Mpc)*GravConst*DiskGravityStellarDiskMass*SolarMass*(DiskGravityStellarDiskScaleHeightR*Mpc+sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)))*fabs(zint)/sqrt(pow(pow(drcyl*LengthUnits,2)+pow((DiskGravityStellarDiskScaleHeightR*Mpc+sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2))),2),3))/sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)));
}

// need to repeat the same functions but with a new r

double func3(double zint)
{
  extern double r2;
  return (-MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)/cosh(fabs(zint)/gScaleHeightz/Mpc)*GravConst*DiskGravityStellarBulgeMass*SolarMass/pow((sqrt(pow(zint,2)+pow(r2,2))+DiskGravityStellarBulgeR*Mpc),2)*fabs(zint)/sqrt(pow(zint,2)+pow(r2,2)));
}

double func4(double zint)
{
  extern double r2;
  return (-MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)/cosh(fabs(zint)/gScaleHeightz/Mpc)*GravConst*DiskGravityStellarDiskMass*SolarMass*(DiskGravityStellarDiskScaleHeightR*Mpc+sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)))*fabs(zint)/sqrt(pow(pow(r2,2)+pow(DiskGravityStellarDiskScaleHeightR*Mpc+sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)),2),3))/sqrt(pow(zint,2)+pow(DiskGravityStellarDiskScaleHeightz*Mpc,2)));
}

// Will be called by qromb to find the pressure at every point in disk.

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n)
{
	static double s;
	static int it;
	int j;
	double del, sum, tnm, x;

	if (n == 1){
		it = 1;
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	}

	tnm = it;
	del = (b-a)/tnm;
	x = a+0.5*del;
	for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
	it *= 2;
	s = 0.5*(s+(b-a)*sum/tnm);
	return s;
} // end trapezoid

#define K 7  // FIXME
FLOAT polint_c[K+1];
FLOAT polint_d[K+1];

/* also called by qromb */
void polint(double xa[],double ya[],int n,double x,double *y,double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	void nrerror(char *);

	dif=fabs(x-xa[1]);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		polint_c[i]=ya[i];
		polint_d[i]=ya[i];
	} // end i for
	
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=polint_c[i+1]-polint_d[i];
			if ( (den=ho-hp) == 0.0 ) fprintf(stderr,"Error in routine POLINT\n");
			den = w/den;
			polint_d[i]=hp*den;
			polint_c[i]=ho*den;
		} // end i for
		*dy=(2*ns < (n-m) ? polint_c[ns+1] : polint_d[ns--]);
		*y += (*dy);
	} // end m for
} // end polint

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP JMAX+1

/* Main integration routine called by DiskPotentialCircularVelocity to find Pressure */
double qromb(double (*func)(double), double a, double b)
{
	if( a == b ) return 0.0;
	double ss,dss,trapzd(double (*func)(double), double a, double b, int n);
  int j;
  double h[JMAXP+1], s[JMAXP+1];
  void polint(double xa[],double ya[],int n,double x,double *y,double *dy),nrerror(char *);

  h[1] = 1.0;
  for (j=1;j<=JMAX;j++){
    s[j] = trapzd(func,a,b,j);
    if( isnan(s[j]) ) ENZO_FAIL("NaN's during pressure integration in GalaxySimulationInitialize");
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j]; 
  }
	/* Print bug report and exit */
  fprintf(stderr,"Too many steps in routine QROMB\n");
  fprintf(stderr,"\t>> drcyl = %"FSYM", z = %"FSYM", z_icm = %"FSYM"\n", drcyl*LengthUnits/Mpc, a/Mpc, b/Mpc);
	fprintf(stderr,"\t>> ss = %"FSYM", dss = %"FSYM"\n", ss, dss);
	ENZO_FAIL("FAILED IN QROMB IN GRID_GALAXYSIMULATIONINIALIZE\n");
	return -1.0;
}
