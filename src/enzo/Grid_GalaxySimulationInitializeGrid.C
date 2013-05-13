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
            float *VelocityUnits, double *MassUnits, FLOAT Time);

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
double DiskPotentialCircularVelocity(FLOAT cellwidth,FLOAT z,FLOAT xpos,FLOAT ypos,FLOAT zpos,FLOAT density,FLOAT &temperature);
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

static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits;

double gScaleHeightR, gScaleHeightz, densicm, MgasScale, Picm;

int grid::GalaxySimulationInitializeGrid(FLOAT DiskRadius,
					 float GalaxyMass,
					 float GasMass,
					 FLOAT DiskPosition[MAX_DIMENSION], 
					 FLOAT ScaleHeightz,
					 FLOAT ScaleHeightR, 
					 float DMConcentration,
					 float DiskTemperature,
					 float InitialTemperature,
					 float UniformDensity,
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

  /* global-scope variables for disk potential functions (FIXME) */
  gScaleHeightR = ScaleHeightR;
  gScaleHeightz = ScaleHeightz;
  densicm = UniformDensity;
  MgasScale = GasMass;
  Picm = kboltz*UniformDensity*InitialTemperature/(0.6*mh);

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
 }

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

 /* Loop over the mesh. */

 float density, dens1, Velocity[MAX_DIMENSION];
 FLOAT temperature, temp1;
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

	density = UniformDensity;
	temperature = temp1 = InitialTemperature;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	/* Find distance from center. */

	r = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		 POW(fabs(y-DiskPosition[1]), 2) +
		 POW(fabs(z-DiskPosition[2]), 2) );
	r = max(r, 0.1*CellWidth[0][0]);

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

	    DiskDensity = (GasMass*SolarMass/(8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))/DensityUnits;   //Code units (rho_0)

	    if( PointSourceGravity > 0 )
	      DiskVelocityMag = gasvel(drad, DiskDensity, ExpansionFactor, GalaxyMass, ScaleHeightR, ScaleHeightz, DMConcentration, Time);
	    else if( DiskGravity > 0 ){
	      CellMass = gauss_mass(drad*LengthUnits,zheight*LengthUnits, xpos*LengthUnits, ypos*LengthUnits, zpos*LengthUnits, inv,
	        DiskDensity*DensityUnits,ScaleHeightR*Mpc, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits);

	      dens1 = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;

	      DiskVelocityMag = DiskPotentialCircularVelocity(CellWidth[0][0], zheight*LengthUnits,xhat[0]*drad*LengthUnits,
	        xhat[1]*drad*LengthUnits,xhat[2]*drad*LengthUnits, dens1, temp1);
			}
      if( PointSourceGravity*DiskGravity != FALSE ) 
	      ENZO_FAIL("Cannot activate both PointSource and Disk gravity options for Isolated Galaxy");

	    if (dim == 0)
	      {
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
	  if (dens1 > density) {
	    density = dens1;
	    if (temp1 == InitialTemperature)
	      temp1 = DiskTemperature;
	    temperature = temp1;
	  }

	} // end: if (r < DiskRadius)
	
	/* Set density. */

	BaryonField[0][n] = density;
	
	if (UseMetallicityField)
	  for (i = 0; i < size; i++)
	    BaryonField[MetalNum][i] = 1.0e-10;
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
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
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

  for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      for (j=0;j<5;j++)
	{
	  yResult[j] = 0.0;
	  for (k=0;k<5;k++)
	    {
	      rot_to_disk(xpos+EvaluationPoints[i]*cellwidth/2.0,ypos+EvaluationPoints[j]*cellwidth/2.0,zpos+EvaluationPoints[k]*cellwidth/2.0,xrot,yrot,zrot,inv);
	      yResult[j] += cellwidth/2.0*Weights[k]*PEXP(-sqrt(POW(xrot,2)+POW(yrot,2))/ScaleHeightR)/POW(cosh(zrot/(2.0*ScaleHeightz)),2);
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

/* 
 *	DISK POTENTIAL CIRCULAR VELOCITY
 */
float DiskPotentialCircularVelocity(FLOAT cellwidth, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT density, FLOAT &temperature)
{

	extern double drcyl;
	double func1(double zint);       //(density times Stellar bulge force)
	double func2(double zint);     //(density times stellar disk force)
	double func3(double zint);       //func1 but for r2
	double func4(double zint);      //func2 but for r2

	double Pressure,Pressure2,zicm,zicm2,zicmf=0.0,zsmall=0.0,
		zicm2f=0.0,zint,FdPdR,FtotR,denuse,rsph,vrot;

	r2=(drcyl+0.01*cellwidth)*LengthUnits;
	rsph=sqrt(pow(drcyl*LengthUnits,2)+pow(z,2));

	/*	Determine zicm: the height above the disk where rho -> rho_ICM,
	 *	use this to find P_icm and dP_icm  */
	if (fabs(drcyl*LengthUnits/Mpc) <= 0.02) {
		zicm=densicm/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc));
		zicm=log(1.0/zicm+sqrt((1.0/pow(zicm,2))-1.0));
		zicm=fabs(zicm*gScaleHeightz*Mpc);


	//	printf("zicm = %g, drcyl = %g\n", zicm/Mpc, drcyl*LengthUnits/Mpc);
		zicm2=densicm/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc));
		zicm2=log(1.0/zicm2+sqrt((1.0/pow(zicm2,2))-1.0));
		zicm2=fabs(zicm2*gScaleHeightz*Mpc);

		Pressure= qromb(func1, fabs(zicm), fabs(z)) + qromb(func2, fabs(zicm), fabs(z));
		Pressure2= qromb(func3, fabs(zicm2), fabs(z)) + qromb(func4, fabs(zicm2), fabs(z));
	}
  else {
    if (fabs(drcyl*LengthUnits/Mpc) <= 0.026) {

			zicmf=densicm/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25
			      /cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)*(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc)))));
			zicm=log(1.0/zicmf+sqrt((1.0/pow(zicmf,2))-1.0));
			zicm=fabs(zicm*gScaleHeightz*Mpc);			
			if (zicmf > 1.0) zicm = 0.0;

			zicm2f=densicm/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25
			      /cosh(r2/gScaleHeightR/Mpc)*(0.5*(1.0+cos(pi*(r2-0.02*Mpc)/(0.006*Mpc)))));
			zicm2=log(1.0/zicm2f+sqrt((1.0/pow(zicm2f,2))-1.0));
			zicm2=fabs(zicm2*gScaleHeightz*Mpc);
			if (zicm2f > 1.0) zicm2 = 0.0;

			if (densicm >= (MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25
			                /cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(z)/gScaleHeightz/Mpc)
			                *(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc))))) 
			    && fabs(z) < zicm) {
				printf("small density zicm = %g, z = %g\n", zicm/Mpc, z/Mpc);
			} // end small density if

			if (fabs(z) < fabs(zicm)) {

				Pressure  = (qromb(func1, fabs(zicm), fabs(z)) + qromb(func2, fabs(zicm), fabs(z)))
				            *(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc))));
    		Pressure2 = (qromb(func3, fabs(zicm2), fabs(z)) + qromb(func4, fabs(zicm2), fabs(z)))
				            *(0.5*(1.0+cos(pi*(r2-0.02*Mpc)/(0.006*Mpc))));
			} // end |z| > |zicm| if

  	} // end r_cyle < .026 if
	} // end r_cyl < .02 if/else

	denuse = density*DensityUnits; // av_den(cellwidth*LengthUnits, z, xpos, ypos, zpos)*DensityUnits; FIXME
	if (Pressure < 0.0 && fabs(drcyl)*LengthUnits/Mpc <= 0.026 && fabs(z) <= fabs(zicm)) {
		fprintf(stderr,"neg pressure:  P = %"FSYM", z = %"FSYM", r = %"FSYM"\n", Pressure, z/Mpc, drcyl*LengthUnits/Mpc);
	}
	if (fabs(drcyl)*LengthUnits/Mpc >= 0.026 || fabs(zicm) <= fabs(z)){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm;
	}
	if (Pressure2 <= 0.0 && Pressure <= 0.0){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm;
	}
	if (Pressure <= 0.0) {
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm;
	}
	if (denuse < densicm) {
		fprintf(stderr,"denuse small:  %"FSYM"\n", denuse);
	}
	temperature=0.6*mh*(Picm+Pressure)/(kboltz*denuse);

	fprintf(stderr,"Temperature, Picm,Pressure,denuse, z,zicm= %"GSYM", %"GSYM", %"GSYM", %"GSYM", %"GSYM", %"GSYM"\n",
		temperature,Picm,Pressure,denuse,z,zicm); // FIXME

	/* Calculate pressure gradient */
	FdPdR = (Pressure2 - Pressure)/(r2-drcyl*LengthUnits)/density; //(av_den(cellwidth*LengthUnits, z, xpos, ypos, zpos)*DensityUnits);

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


double *vector(int nl,int nh)
{
  double *v;
  v = new double[(nh-nl+1)];
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}

void free_vector(double *v,int nl,int nh)
{
  delete [] &v[nl];
}

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

/* also called by qromb */
void polint(double xa[],double ya[],int n,double x,double *y,double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d,*vector(int,int);
	void nrerror(char *),free_vector(double *, int, int);

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	} // end i for
	
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0 ) fprintf(stderr,"Error in routine POLINT\n");
			den = w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		} // end i for
		*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]);
		*y += (*dy);
	} // end m for

	free_vector(d,1,n);
	free_vector(c,1,n);
} // end polint

#define EPS 1.0e-6
#define JMAX 40
#define JMAXP JMAX+1
#define K 7

/* Main integration routine called by DiskPotentialCircularVelocity to find Pressure */
double qromb(double (*func)(double), double a, double b)
{
	double ss,dss,trapzd(double (*func)(double), double a, double b, int n);
  int j;
  double h[JMAXP+1], s[JMAXP+1];
  void polint(double xa[],double ya[],int n,double x,double *y,double *dy),nrerror(char *);

  h[1] = 1.0;
  for (j=1;j<=JMAX;j++){
    s[j] = trapzd(func,a,b,j);
    if( s[j] != s[j] ) ENZO_FAIL("NaN's during pressure integration in GalaxySimulationInitialize");
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j]; 
  }
  fprintf(stderr,"Too many steps in routine QROMB\n");
  fprintf(stderr,"\t>> drcyl= %"FSYM", z = %"FSYM"\n", drcyl*LengthUnits/Mpc, a/Mpc); 
	return -1.0;
}
