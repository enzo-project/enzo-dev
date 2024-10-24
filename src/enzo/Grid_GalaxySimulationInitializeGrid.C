/***********************************************************************
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
#include <assert.h>
#include "preincludes.h" 
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "phys_constants.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define kboltzKeV (8.617e-8)    // keV per K
#define mu (0.6)
#define CM_PER_KM (1.0e5)
#define CM_PER_KPC (3.0856e21)
#define KEV_PER_ERG (6.242e8)

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                      float *TemperatureUnits, float *TimeUnits,
                      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* Internal routines */
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float&  HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest);
double gasvel(FLOAT radius, double DiskDensity, FLOAT ExpansionFactor, 
             double GalaxyMass, double ScaleHeightR, double ScaleHeightz, 
             double DMConcentration, FLOAT Time);
double gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, 
                 FLOAT inv [3][3], double DiskDensity,
                 double ScaleHeightR, double ScaleHeightz,
                 FLOAT cellwidth);
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot,
                 FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);
double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2);

double DarkMatterMassEnclosed(FLOAT R);

/* struct to carry around data required for circumgalactic media
   if we need to generate radial profiles of halo quantities via 
   numerical integration */
struct CGMdata {
  double *n_rad, *T_rad, *rad, *press;
  int nbins;
  double R_inner, R_outer, dr;

  CGMdata(int n_bins) {
    
    nbins = n_bins;
  
    n_rad = new double[nbins];
    T_rad = new double[nbins];
    rad   = new double[nbins];

    for(int i=0; i<nbins; i++) n_rad[i] = T_rad[i] = rad[i] = -1.0;
  }
  
  ~CGMdata() {
    if (n_rad) delete[] n_rad;
    if (T_rad) delete[] T_rad;
    if (rad) delete[] rad;
  }
};

/* Internal Routines for CGM setup */
double HaloGasDensity(FLOAT R, struct CGMdata&);
double HaloGasTemperature(FLOAT R, struct CGMdata&);

/* Internal Routines for DiskGravity Setup */
double DiskGravityCircularVelocity(FLOAT rsph, FLOAT rcyl, FLOAT z);
double DiskGravityStellarAccel(FLOAT rcyl, FLOAT z);
double DiskGravityBulgeAccel(FLOAT rsph);

/* Stuff for the Tonnesen & Bryan 09 ressurection */
void DiskForceBalance(FLOAT cellwidth, FLOAT z, double density, struct CGMdata& CGM_data,
                          double &temperature, double &rotvel);
double trapzd(double (func)(), double a, double b, int n);
double qromb(double (*func)(double), double a, double b);
void polint(double xa[],double ya[],int n,double x,double *y,double *dy);
static FLOAT r2;

static float DensityUnits, LengthUnits, TemperatureUnits = 1,
             TimeUnits, VelocityUnits, MassUnits;

double gScaleHeightR, gScaleHeightz, densicm, MgasScale, Picm,
       TruncRadius, SmoothRadius, SmoothLength,Ticm;

/* Global variables (within this file) for circumgalactic medium setup 
   (also used a bit for disk potential setup) */
int GalaxySimulationGasHalo, EquilibrateChem;
double GalaxySimulationGasHaloScaleRadius,
  GalaxySimulationGasHaloDensity, GalaxySimulationGasHaloDensity2,
  GalaxySimulationGasHaloTemperature, GalaxySimulationGasHaloAlpha,
  GalaxySimulationGasHaloZeta, GalaxySimulationGasHaloZeta2,
  GalaxySimulationGasHaloCoreEntropy, 
  GalaxySimulationGalaxyMass, GalaxySimulationDMConcentration,
  GalaxySimulationGasHaloMetallicity,
  GalaxySimulationDiskMetallicityEnhancementFactor,
  GalaxySimulationGasHaloRatio;

/* declarations for a bunch of functions needed to generate radial profiles
   of halo quantities via numerical integration - see the actual functions for
   descriptions. */
double sigmoid(FLOAT r, FLOAT r0, double k, double y0, double y_off);
double halo_S_of_r(FLOAT r);
double halo_S_of_r(FLOAT r, grid* Grid);
double halo_dSdr(FLOAT r, double n);
double halo_dn_dr(FLOAT r, double n);
double halo_dP_dr(FLOAT r, double P, grid* Grid);
double sigmoid_halo_dP_dr(FLOAT r, double P, double log_r0,
			  double k, double y0, double y_off);
double halo_g_of_r(FLOAT r);
double halo_mod_g_of_r(FLOAT r);
double halo_mod_DMmass_at_r(FLOAT r);
void halo_init(struct CGMdata& CGM_data, grid* Grid, FLOAT Rstop=-1, int GasHalo_override=0);

static FLOAT rcyl;

int grid::GalaxySimulationInitializeGrid(double DiskRadius,
           double GalaxyMass,
           double GasMass,
           FLOAT DiskPosition[MAX_DIMENSION], 
           double ScaleHeightz,
           double ScaleHeightR,
           double GalaxyTruncationRadius,
           double DiskDensityCap, 
           double DMConcentration,
           double DiskTemperature,
           int    DiskPressureBalance,
           double InitialTemperature,
           double UniformDensity,
           int   EquilChem,
	   int   GasHalo,
           double GasHaloScaleRadius,
           double GasHaloDensity,
           double GasHaloDensity2,
           double GasHaloTemperature,
           double GasHaloAlpha,
           double GasHaloZeta,
           double GasHaloZeta2,
           double GasHaloCoreEntropy,
	         double GasHaloRatio,           
           double GasHaloMetallicity,
           int   UseHaloRotation,
           double RotationScaleVelocity,
           double RotationScaleRadius,
           double RotationPowerLawIndex,
           double DiskMetallicityEnhancementFactor,
           double AngularMomentum[MAX_DIMENSION],
           double UniformVelocity[MAX_DIMENSION], 
           int UseMetallicityField, 
           FLOAT GalaxySimulationInflowTime,
           double GalaxySimulationInflowDensity,
           int level,
           double GalaxySimulationCR
          )
{
  /* declarations */

  int dim, i, j, k, m, field, disk, size, MetalNum, MetalIaNum, vel;
  int MetalIINum, MetalAGBNum, MetalNSMNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
    H2IINum, DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
  double DiskDensity, DiskVelocityMag;
  int CRNum, DensNum;

  /* global-scope variables for disk potential functions (would be better if not global) */

  gScaleHeightR = ScaleHeightR;
  gScaleHeightz = ScaleHeightz;
  densicm = UniformDensity;  // gas density if no halo is used
  MgasScale = GasMass;
  Ticm = InitialTemperature;  // gas temperature if no halo is used
  Picm = kboltz*UniformDensity*Ticm/(mu*mh);  // gas pressure if no halo is used
  TruncRadius = GalaxyTruncationRadius;
  SmoothRadius = TruncRadius*.02/.026;
  SmoothLength = TruncRadius - SmoothRadius;

  /* set all of the gas halo quantities to variables that are global
     within this file, so the calls to HaloGasDensity() and HaloGasTemperature()
     are unchanged. */
  EquilibrateChem = EquilChem;
  GalaxySimulationGasHalo = GasHalo;   // integer, >= 0
  GalaxySimulationGasHaloScaleRadius = GasHaloScaleRadius;  // in mpc
  GalaxySimulationGasHaloDensity = GasHaloDensity; // in grams/cm^3
  GalaxySimulationGasHaloDensity2 = GasHaloDensity2;
  GalaxySimulationGasHaloTemperature = GasHaloTemperature;  // in Kelvin
  GalaxySimulationGasHaloAlpha = GasHaloAlpha;  // power-law index; unitless
  GalaxySimulationGasHaloZeta = GasHaloZeta;
  GalaxySimulationGasHaloZeta2 = GasHaloZeta2;
  GalaxySimulationGasHaloCoreEntropy = GasHaloCoreEntropy;  // power-law index; unitless
  GalaxySimulationGasHaloRatio = GasHaloRatio; // ratio of cooling time to freefall time
  GalaxySimulationGalaxyMass = GalaxyMass;
  GalaxySimulationDMConcentration = DMConcentration;
  GalaxySimulationGasHaloMetallicity = GasHaloMetallicity; // Zsun
  GalaxySimulationDiskMetallicityEnhancementFactor = DiskMetallicityEnhancementFactor; // w.r.t to halo
 
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
  if (UseMHD) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
  }
  if(HydroMethod == MHD_RK ){
    FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
  }
  if (UsePoissonDivergenceCleaning) {
    FieldType[NumberOfBaryonFields++] = Phi_pField;
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
  if (StarMakerTypeIaSNe || StarFeedbackTrackMetalSources)
    FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;
  if (StarFeedbackTrackMetalSources) {
    FieldType[MetalIINum = NumberOfBaryonFields++] = MetalSNIIDensity;
    FieldType[MetalAGBNum = NumberOfBaryonFields++] = MetalAGBDensity;
    FieldType[MetalNSMNum = NumberOfBaryonFields++] = MetalNSMDensity;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;

  /* Set various units. */

  double CriticalDensity = 1, BoxLength = 1;
  FLOAT a, dadt, ExpansionFactor = 1;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);
    CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc_cm^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc_cm
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

  /* Scale gas halo rotation quantities to code units.
   * gas halo rotation variable are NOT global */
  RotationScaleVelocity *= CM_PER_KM; // km/s to cm/s
  RotationScaleVelocity /= LengthUnits/TimeUnits; // cm/s to code length/code time
  RotationScaleRadius *= CM_PER_KPC;  // kpc to cm
  RotationScaleRadius /= LengthUnits;  // cm to code length

  /*  initializes halo radius, density, temperature profiles 
      for circumgalactic medium if needed (i.e., for CGM profiles that
      require integration to get quantities we care about. 
      Assumes disk is located at the center of the domain. */

  FLOAT far_left, far_right, largest_rad;
  
  far_left = DomainLeftEdge[0];
  far_right = DomainRightEdge[0];
  
  for (int i=1; i<GridRank; ++i) {
    if (DomainLeftEdge[i] < far_left)
      far_left = DomainLeftEdge[i];
    if (DomainRightEdge[i] > far_right)
      far_right = DomainRightEdge[i];
  }

  largest_rad = sqrt(3) * (far_right - far_left) / 2.0 * LengthUnits;

  struct CGMdata CGM_data(8192);
  halo_init(CGM_data, this, largest_rad);
  if (debug) printf("Made halo profile\n");

  // for (int i=0; i<CGM_data.nbins; ++i)
  //   printf("%g %g %g %g\n", CGM_data.rad[i], CGM_data.n_rad[i], CGM_data.T_rad[i], CGM_data.press[i]);
  
  /* compute size of fields */
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */
  this->AllocateGrids();

  /* I'm commenting this out because the metal field should
     be set during grid initialization rather than just setting
     it as a constant color field. -- DWS */
  // /* set metals to small value */
  //  if (UseMetallicityField)
  //    for (i = 0; i < size; i++)
  //      BaryonField[MetalNum][i] = 1.0e-10;
 
  /* Loop over the mesh. */
  double density, disk_dens;
  double halo_vmag, disk_vel[MAX_DIMENSION], Velocity[MAX_DIMENSION];
  double temperature, disk_temp, init_temp, initial_metallicity;
  FLOAT r_sph, x, y = 0, z = 0;
  int n = 0, iter;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

	if (UseMetallicityField) {
	  /* Set a background metallicity value that will scale with density.
	     If the cell is in the disk, this wifll be increased by a factor
	     of 3.  This should really be a parameter that is read in -- DWS */ 
	  initial_metallicity = GalaxySimulationGasHaloMetallicity;
	}

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	for (dim = 0; dim < MAX_DIMENSION; dim++){
	  Velocity[dim] = 0;
	  disk_vel[dim] = 0;
  }

	/* Find distance from center. */

	r_sph = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		     POW(fabs(y-DiskPosition[1]), 2) +
		     POW(fabs(z-DiskPosition[2]), 2) );
	r_sph = max(r_sph, 0.1*CellWidth[0][0]);
    
	/*
	  r_cyl = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
	  POW(fabs(y-DiskPosition[1]), 2) );
	*/

	density = HaloGasDensity(r_sph, CGM_data)/DensityUnits;
	temperature = init_temp = HaloGasTemperature(r_sph, CGM_data);
  disk_temp = DiskTemperature;
	
	FLOAT xpos, ypos, zpos, rsph, zheight, theta; 
	double CellMass;
	FLOAT rp_hat[3];
	FLOAT yhat[3];

	/* Loop over dims if using Zeus (since vel's face-centered). */

	for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
	     dim++) {

	  /* Compute position. */

	  xpos = x-DiskPosition[0]-(dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	  ypos = y-DiskPosition[1]-(dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	  zpos = z-DiskPosition[2]-(dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

	  /* Compute z and r_perp (AngularMomentum is angular momentum 
	     and must have unit length). */    

	  /* magnitude of z = r.L in L direction */

	  zheight = AngularMomentum[0]*xpos + 
	    AngularMomentum[1]*ypos +
	    AngularMomentum[2]*zpos;

	  /* position in plane of disk */

	  rp_hat[0] = xpos - zheight*AngularMomentum[0];
	  rp_hat[1] = ypos - zheight*AngularMomentum[1];
	  rp_hat[2] = zpos - zheight*AngularMomentum[2];
	  rcyl = sqrt(rp_hat[0]*rp_hat[0] + rp_hat[1]*rp_hat[1] + rp_hat[2]*rp_hat[2]);

	  /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	  rp_hat[0] = rp_hat[0]/rcyl;
	  rp_hat[1] = rp_hat[1]/rcyl;
	  rp_hat[2] = rp_hat[2]/rcyl;
      

	  /* If requested, calculate velocity for CGM halo.
	   * Will be replaced wtih disk velocity later if appropriate */
	  if (UseHaloRotation){
	    /* polar angle as measured from the angular momentum vector*/
	    theta = acos(zheight/r_sph);

	    halo_vmag = RotationScaleVelocity // code units
	      * POW(r_sph/RotationScaleRadius, 
		    RotationPowerLawIndex);

	    if (r_sph <= RotationScaleRadius)
	      halo_vmag = RotationScaleVelocity;

	    halo_vmag *= sin(theta)*sin(theta);
	  
	    /* Cylindrical velocity */
	    Velocity[0] = halo_vmag * (AngularMomentum[1]*rp_hat[2] -
				       AngularMomentum[2]*rp_hat[1]);
	    Velocity[1] = halo_vmag * (AngularMomentum[2]*rp_hat[0] -
				       AngularMomentum[0]*rp_hat[2]);
	    Velocity[2] = halo_vmag * (AngularMomentum[0]*rp_hat[1] -
				       AngularMomentum[1]*rp_hat[0]);
	  }

	  disk_dens = 0.0;
	  if (r_sph < DiskRadius) {

	    /* Beyond truncation radius */
	    if( fabs(rcyl*LengthUnits/Mpc_cm) > TruncRadius ){
	      disk_dens = 0.0;
	      break;
	    }

	    /* Find another vector perpendicular to r_perp and AngularMomentum */

	    yhat[0] = AngularMomentum[1]*rp_hat[2] - AngularMomentum[2]*rp_hat[1];
	    yhat[1] = AngularMomentum[2]*rp_hat[0] - AngularMomentum[0]*rp_hat[2];
	    yhat[2] = AngularMomentum[0]*rp_hat[1] - AngularMomentum[1]*rp_hat[0];

	    /* generate rotation matrix */
	    FLOAT inv[3][3],temp;
	    int i,j;

	    // matrix of basis vectors in coordinate system defined by the galaxy
	    inv[0][0] = rp_hat[0];
	    inv[0][1] = yhat[0];
	    inv[0][2] = AngularMomentum[0];
        
	    inv[1][0] = rp_hat[1];
	    inv[1][1] = yhat[1];
	    inv[1][2] = AngularMomentum[1];
        
	    inv[2][0] = rp_hat[2];
	    inv[2][1] = yhat[2];
	    inv[2][2] = AngularMomentum[2];

	    // Matrix is orthogonal by construction so inverse = transpose
	    for (i=0;i<3;i++)
	      for (j=i+1;j<3;j++){
          temp = inv[i][j];
          inv[i][j] = inv[j][i];
          inv[j][i] = temp;
	      }

	    DiskDensity = (GasMass * SolarMass
			   / (8.0*pi*ScaleHeightz*Mpc_cm*POW(ScaleHeightR*Mpc_cm,2.0)))
	      / DensityUnits;   //Code units (rho_0) 

	    CellMass = gauss_mass(rcyl*LengthUnits, zheight*LengthUnits,
				  xpos*LengthUnits, ypos*LengthUnits,
				  zpos*LengthUnits, inv, 
				  DiskDensity*DensityUnits,
				  ScaleHeightR*Mpc_cm, ScaleHeightz*Mpc_cm, 
				  CellWidth[0][0]*LengthUnits);

	    disk_dens = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;

	    if ((disk_dens > DiskDensityCap) && (DiskDensityCap > 0))
	      disk_dens = DiskDensityCap;

	    /* Inside CGM */
	    if (disk_dens < density)
	      break;

	    //
	    // calculate velocity
	    //

	    if (PointSourceGravity > 0 )
	      DiskVelocityMag = gasvel(rcyl, DiskDensity, ExpansionFactor,
				       GalaxyMass, ScaleHeightR,
				       ScaleHeightz, DMConcentration, Time);

	    else if( DiskGravity > 0 ) {

        if (!DiskPressureBalance)
          DiskVelocityMag = DiskGravityCircularVelocity(r_sph*LengthUnits,
                                                        rcyl*LengthUnits,
                                                        zheight*LengthUnits)
                            /VelocityUnits;
        else // will also set disk_temp to be not isothermal and rotation will account for this
          DiskForceBalance(CellWidth[0][0], zheight*LengthUnits, disk_dens, CGM_data, disk_temp, DiskVelocityMag);
      }
        
	    if (PointSourceGravity*DiskGravity != FALSE ) 
	      ENZO_FAIL("Cannot activate both PointSource and Disk gravity options for Isolated Galaxy");

	    /* Compute velocty: L x r_perp. */
	    if (dim == 0 || dim == 1)
	      disk_vel[0] = DiskVelocityMag*(AngularMomentum[1]*rp_hat[2] -
					     AngularMomentum[2]*rp_hat[1]);
	    if (dim == 0 || dim == 2)
	      disk_vel[1] = DiskVelocityMag*(AngularMomentum[2]*rp_hat[0] -
					     AngularMomentum[0]*rp_hat[2]);
	    if (dim == 0 || dim == 3)
	      disk_vel[2] = DiskVelocityMag*(AngularMomentum[0]*rp_hat[1] -
					     AngularMomentum[1]*rp_hat[0]);

	  } // end: if (r_sph < DiskRadius)

	  /* Replace CGM ("Halo") defaults with disk if dense enough; i.e.
	   * replace 'density', 'temperature', 'initial_metallicity', and
	   * 'Velocity' (which are currently set to CGM values) with their
	   * appropriate disk values */
       
	  if (disk_dens > density && fabs(rcyl*LengthUnits/Mpc_cm) <= TruncRadius){
        
	    density = disk_dens;
	    temperature = disk_temp;
        
	    /* Here we're setting the disk to be X times more enriched -- DWS */
	    if( UseMetallicityField )
	      initial_metallicity *= GalaxySimulationDiskMetallicityEnhancementFactor;
          
	    /* Replace default/CGM velocity with disk velocity */
	    Velocity[0] = disk_vel[0];
	    Velocity[1] = disk_vel[1];
	    Velocity[2] = disk_vel[2];
	  }

	} // end: loop over dims 

	/* Set density. */

	BaryonField[0][n] = density;

	if (UseMetallicityField) {
	  BaryonField[MetalNum][n] = initial_metallicity 
	    * CoolData.SolarMetalFractionByMass 
	    * density;
	}

	/* This should probably be scaled with density in some way to be
	   a proper metallicity -- DWS (loop redundancy addressed by CEK) */
	if (StarMakerTypeIaSNe || StarFeedbackTrackMetalSources)
	  BaryonField[MetalIaNum][n] = tiny_number;

  if (StarFeedbackTrackMetalSources){
    BaryonField[MetalIINum][n] = tiny_number;
    BaryonField[MetalAGBNum][n] = tiny_number;
    BaryonField[MetalNSMNum][n] = tiny_number;
  }
   
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits / ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];

	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);

	if (BaryonField[1][n] <= 0.0)
	  printf("G_GSIC: negative or zero energy  n = %"ISYM"  temp = %"FSYM"   e = %"FSYM"\n",
		 n, temperature, BaryonField[1][n]);


	// Set multispecies fields!
	// this attempts to set them such that species conservation is maintained,
	// using the method in CosmologySimulationInitializeGrid.C
	if(MultiSpecies){
	  if (MultiSpecies == 3)
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		       BaryonField[DINum][n], BaryonField[DIINum][n], BaryonField[HDINum][n]);
	  else if (MultiSpecies == 2) {
	    float temp;
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		       temp, temp, temp);
	  }
	  else {
	    float temp;
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       temp, temp, temp,
		       temp, temp, temp);
	  }
	} // if(MultiSpecies)
	
      } // end loop over grids
    
  if( CRModel )
    BaryonField[CRNum][n] = BaryonField[DensNum][n] * GalaxySimulationCR;

  return SUCCESS;

} // end Grid::GalaxySimulationInitializeGrid

/*
* Initialization routines
* Order:
*   halo mass (DarkMatterMassEnclosed)
*   cell mass (gauss_mass, rot_to_disk)
*   disk velocity (gas_vel+DiskGravity OR DiskPotential+findZicm+DiskForceBalance etc)
*   chemistry (setup_chem, bilinear_interp)
*   CGM profile
*/


/* halo galaxy mass at a given radius, using user-defined global parameters for galaxy
   quantities and assuming that all halo mass is in an NFW halo.  This is not totally
   correct near the center of the halo, but since we're using it for the CGM initialization 
   and are dealing with radii that aren't particularly near the center of the halo, this 
   approximation is probably fine. 

   Input is the radius in CGS units; output is the enclosed mass at that radius in CGS units.
*/
double DarkMatterMassEnclosed(FLOAT r){

  double M_within_r;

  if (DiskGravityDarkMatterUseNFW) {
    double M, C, R200, rho_0, Rs;
    double rho_crit = 1.8788e-29*0.49;
    
    // GSDarkMatterConcentration is the same as DiskGravityDarkMatterConcentration
    // if the latter is in use, same with GSGalaxyMass & DGDarkMatterMass

    M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
    C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
    
    R200 = POW(3.0/(4.0*pi)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
    Rs = R200/C;  // scale radius of NFW halo in CGS
    rho_0 = 200.0*POW(C,3)/3.0/(log(1.0+C) - C/(1.0+C))*rho_crit;  // rho_0 for NFW halo in CGS

    // mass w/in radius R
    M_within_r = 4.0*pi*rho_0*POW(Rs,3.0)*(log((Rs+r)/Rs) - r/(Rs+r));

  } else if (DiskGravityDarkMatterUseB95) {
    FLOAT R0 = DiskGravityDarkMatterR*Mpc_cm;
    FLOAT x = r/R0*LengthUnits;
	  double M0 = pi*DiskGravityDarkMatterDensity*R0*R0*R0;

	  M_within_r = M0*(-2.0*atan(x)+2.0*log(1+x)+log(1.0+x*x));

  } else {
    ENZO_FAIL("This GalaxySimulationGasHalo option requires a DiskGravityDarkMatter profile.")
  }

  return M_within_r;
  
}

// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature
double gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], 
                  double DiskDensity, double ScaleHeightR, double ScaleHeightz, FLOAT cellwidth)
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
                    if( rrot/Mpc_cm < SmoothRadius )
                        yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz);
                    else if( rrot/Mpc_cm < TruncRadius )
                        yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz)
                                        *0.5*(1.0+cos(pi*(rrot-SmoothRadius*Mpc_cm)/(SmoothLength*Mpc_cm)));
                } // end disk gravity if

            }
            xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
        }
        Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }  
    Mass *= DiskDensity;
    return Mass;
}

//Finds coordinates in rotated coordinate system; used by gauss_mass
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3])
{
  xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
  yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
  zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];
}



//
// Disk velocity with PointSourceGravity
//
double gasvel(FLOAT radius, double DiskDensity, FLOAT ExpansionFactor, 
              double GalaxyMass, double ScaleHeightR, double ScaleHeightz,
              double DMConcentration, FLOAT Time)
{

 double OMEGA=OmegaLambdaNow+OmegaMatterNow;                 //Flat Universe

 double r = radius*LengthUnits/100;    // Radius [m]

 double M_200 = GalaxyMass*SolarMass/1000.0;      // Virial Mass [kg]

 double H = sqrt(HubbleConstantNow*100*HubbleConstantNow*100*(OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor,-3)-(OMEGA-1.)*POW(ExpansionFactor,-2)));                                

 double r_200 = (1.63e-2*POW(GalaxyMass,1.0/3.0)*POW((OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor, -3)-(OMEGA-1.0)*POW(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*POW(H,-2.0/3.0)*POW(100,2.0/3.0))*Mpc_cm/1.0e5;
 //virial radius [m]: M_200/M_Solar = GalaxyMass

 double M_gas, M_DM, M_Tot, Acc, V_Circ;
 double f_C = log(1.0+DMConcentration)-DMConcentration/(1.0+DMConcentration);
 double r_s = r_200/DMConcentration;  //[m]

 // Mass of gas disk and DM at given radius

  M_gas=8.0*M_PI*ScaleHeightz*Mpc_cm/100*ScaleHeightR*Mpc_cm/100*ScaleHeightR*Mpc_cm/100*DiskDensity*DensityUnits*1000*PEXP(-r/(ScaleHeightR*Mpc_cm/100))*(PEXP(r/(ScaleHeightR*Mpc_cm/100))-r/(ScaleHeightR*Mpc_cm/100)-1.0);

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

  // Force per unit mass on disk (i.e. acceleration) [ms-2]

  Acc=((GravConst/1000.0)*M_Tot)/(r*r);

 // Magnitude of Circular Velocity of disk 

  V_Circ = sqrt(r*Acc)*100;       //cms-1

  return (V_Circ/VelocityUnits);  //code units
}


//
// Disk velocity with DiskGravity. Arguments should be in cgs
//
double DiskGravityStellarAccel(FLOAT rcyl, FLOAT z){
    // return (cylindrical) radial component of acceleration ...? in cgs?
    // Potential from  Miyamoto & Nagai 75

    double accelcylR;

    accelcylR = GravConst * DiskGravityStellarDiskMass*SolarMass * rcyl
              / sqrt(POW(
                    POW(rcyl,2)
                  + POW(DiskGravityStellarDiskScaleHeightR*Mpc_cm
                        + sqrt(POW(z,2)
                             + POW(DiskGravityStellarDiskScaleHeightz*Mpc_cm,2)),
                        2),
                     3));

    return accelcylR;
}

double DiskGravityBulgeAccel(FLOAT rsph) { // cgs arguments
    double accelsph;

    accelsph = GravConst * DiskGravityStellarBulgeMass*SolarMass
             / POW(rsph + DiskGravityStellarBulgeR*Mpc_cm,2);

    return accelsph;
}

double DiskGravityCircularVelocity(FLOAT rsph, FLOAT rcyl, FLOAT z) {
    double acc=0, velmag;
    
    if (DiskGravityDarkMatterUseNFW || DiskGravityDarkMatterUseB95)
      acc = halo_g_of_r(rsph);

    acc += DiskGravityStellarAccel(rcyl, z)
         + DiskGravityBulgeAccel(rsph);

    velmag = sqrt(acc*rcyl); // cgs
    return velmag;
}

/* -------------------- BEGINNING of Functions for thermal pressure balance in the disk --------------------
  Assumes Burkert 95 DM potential. Originally designed for galaxies in an ICM wind. See Tonnesen & Bryan 09.
*/

float DiskPotentialGasDensity(FLOAT r,FLOAT z){
/*
 *  computes gas density within galaxy disk, according to eq
 *
 *    (Mgas/8*pi*a^2*b)*sech(r/a)*sech*(z/b)
 *
 *  Smoothed by a cosine fcn beyond SmoothRadius
 *
 *  Parameteres:
 *  ------------
 *    r - cylindrical radius (code units)
 *    z - cylindrical height (code units)
 *
 *  Returns: density (in grams/cm^3)
 *
 */
  double density = MgasScale*SolarMass/(8.0*pi*POW(gScaleHeightR*Mpc_cm,2)*gScaleHeightz*Mpc_cm);
  density /= (cosh(r*LengthUnits/gScaleHeightR/Mpc_cm)*cosh(z*LengthUnits/gScaleHeightz/Mpc_cm));

  if(fabs(r*LengthUnits/Mpc_cm) > SmoothRadius && fabs(r*LengthUnits/Mpc_cm) <= TruncRadius)
    density *= 0.5*(1.0+cos(pi*(r*LengthUnits-SmoothRadius*Mpc_cm)/(SmoothLength*Mpc_cm)));
  return density;
} // end DiskPotentialGasDensity


double findZicm(FLOAT r, struct CGMdata& CGM_data){
  /*  
   *  Finds the height above the disk plane where the disk gas density
   *  matches the halo's gas density (using bisection)
   *
   *  Parameters:
   *  -----------
   *    r - cylindrical radius (code units)
   *
   *  Returns: zicm, edge of disk, (code units)
   */

  static const double X_TOL = 1e-7*Mpc_cm/LengthUnits; // sub pc resolution
  static const int MAX_ITERS = 50; int iters=0;

  double z_lo = 0.0,z_hi = 0.01*Mpc_cm/LengthUnits,z_new,f_lo,f_hi,f_new;
  f_hi = DiskPotentialGasDensity(r,z_hi) - HaloGasDensity(sqrt(r*r+z_hi*z_hi), CGM_data); // -ve
  f_lo = DiskPotentialGasDensity(r,z_lo) - HaloGasDensity(sqrt(r*r+z_lo*z_lo), CGM_data); // +ve

  if(f_lo < 0.0) return 0.0; // beyond the disk
  if(f_hi > 0.0) ENZO_FAIL("ERROR IN GALAXY INITIALIZE: HALO IS UNDER-PRESSURIZED");

  while(iters++ < MAX_ITERS ){

    z_new = (z_hi+z_lo)/2.0;
    f_new = DiskPotentialGasDensity(r,z_new)
            - HaloGasDensity(sqrt(r*r+z_new*z_new), CGM_data);

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

void DiskForceBalance(FLOAT cellwidth, FLOAT z, double density, struct CGMdata& CGM_data, double &temperature, double &rotvel)
{
  /* 
   *  DISK POTENTIAL CIRCULAR VELOCITY
   *
   *      Sets disk temperature & circular velocity (in code units) given height z
   *      and rcyl (radius in plane) in disk. Velocity includes the effect of
   *      thermal pressure and so is not just the sqrt(G M/rcyl).
   *      Note that for historical reasons rcyl is an external. *
   */

  extern FLOAT rcyl;                    // will be in code units
  double PbulgeComp1(double zint);       // (density times Stellar bulge force)
  double PbulgeComp2(double zint);       // same but for r2 (3D distance)
  double PstellarComp1(double zint);     // (density times stellar disk force)
  double PstellarComp2(double zint);     // same but for r2 (3D distance plane)
  double PDMComp1(double zint);          // (density times dark matter halo force)
  double PDMComp2(double zint);          // same but for r2 (3D distance plane)

  double Pressure,Pressure2,zicm,zicm2,zicmf=0.0,zsmall=0.0,
    zicm2f=0.0,zint,FdPdR,FtotR,denuse,rsph,vrot,bulgeComp,rsph_icm;

  /* Distances in cgs */
  r2 = (rcyl+0.01*cellwidth)*LengthUnits;  // in plane radius
  rsph = sqrt(POW(rcyl*LengthUnits,2)+POW(z,2)); // 3D radius

  /*  Determine zicm: the height above the disk where rho -> rho_ICM,
   *  use this to find P_icm and dP_icm  */

  if (fabs(rcyl*LengthUnits/Mpc_cm) <= SmoothRadius) {

    zicm  = findZicm(rcyl, CGM_data)*LengthUnits;
    zicm2 = findZicm(r2/LengthUnits, CGM_data)*LengthUnits; // r2 = rcyl + delta

    if( fabs(z) < fabs(zicm) ){

      /* Integrate the density times force to get pressure.  Do this
         at two different locations to get a numerical gradient. */

      bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ? 
       0.0 : qromb(PbulgeComp1, fabs(zicm), fabs(z)));
      Pressure  = bulgeComp + qromb(PstellarComp1, fabs(zicm), fabs(z));
      Pressure += qromb(PDMComp1, fabs(zicm), fabs(z));

      bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ? 
       0.0 : qromb(PbulgeComp2, fabs(zicm2), fabs(z)));
      Pressure2  = bulgeComp + qromb(PstellarComp2, fabs(zicm2), fabs(z));
      Pressure2 += qromb(PDMComp2, fabs(zicm2), fabs(z));

    }  // end |z| < |zicm| if

  }  else {

    if (fabs(rcyl*LengthUnits/Mpc_cm) <= TruncRadius ) {

      zicm  = findZicm(rcyl, CGM_data)*LengthUnits;
      zicm2 = findZicm(r2/LengthUnits, CGM_data)*LengthUnits;


      if (fabs(z) < fabs(zicm)) {
        
        bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ?
              0.0 : qromb(PbulgeComp1, fabs(zicm), fabs(z)));
        Pressure  = (bulgeComp + qromb(PDMComp1, fabs(zicm), fabs(z)) + qromb(PstellarComp1, fabs(zicm), fabs(z)))
          *(0.5*(1.0+cos(pi*(rcyl*LengthUnits-SmoothRadius*Mpc_cm)/
            (SmoothLength*Mpc_cm))));

        bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ?
              0.0 : qromb(PbulgeComp2, fabs(zicm2), fabs(z)));
        Pressure2 = (bulgeComp + qromb(PDMComp2, fabs(zicm2), fabs(z)) + qromb(PstellarComp2, fabs(zicm2), fabs(z)))
          *(0.5*(1.0+cos(pi*(r2-SmoothRadius*Mpc_cm)/(SmoothLength*Mpc_cm))));

      } // end |z| < |zicm| if

    } // end r_cyle < TruncRadius if

  } // end r_cyl < SmoothRadius if/else

  denuse = density*DensityUnits; 

  if (Pressure < 0.0 && fabs(rcyl)*LengthUnits/Mpc_cm <= TruncRadius && fabs(z) <= fabs(zicm)) {
    fprintf(stderr,"neg pressure:  P = %"FSYM", z = %"FSYM", r = %"FSYM"\n", Pressure, z/Mpc_cm, rcyl*LengthUnits/Mpc_cm);
  }
  if (fabs(rcyl)*LengthUnits/Mpc_cm >= TruncRadius || fabs(zicm) <= fabs(z)){
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph/LengthUnits, CGM_data);
  }
  if (Pressure2 <= 0.0 && Pressure <= 0.0){
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph/LengthUnits, CGM_data);
  }
  if (Pressure <= 0.0) {
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph/LengthUnits, CGM_data);
  }
  if (denuse < HaloGasDensity(rsph/LengthUnits, CGM_data)) {
    fprintf(stderr,"denuse small:  %"FSYM"\n", denuse);
  }
  
  rsph_icm = sqrt(rcyl*rcyl+POW(zicm/LengthUnits,2));  // code units
  Picm = HaloGasDensity(rsph_icm, CGM_data)*kboltz*HaloGasTemperature(rsph_icm, CGM_data)/(mu*mh);
  temperature=mu*mh*(Picm+Pressure)/(kboltz*denuse);

  /* Calculate pressure gradient */

  FdPdR = (Pressure2 - Pressure)/(r2-rcyl*LengthUnits)/density; 

  /* Calculate Gravity = Fg_DM + Fg_StellarDisk + Fg_StellaDiskGravityStellarBulgeR */
  
  FtotR  = (-pi)*GravConst*DiskGravityDarkMatterDensity*
            POW(DiskGravityDarkMatterR*Mpc_cm,3)/POW(rsph,3)*rcyl*LengthUnits
            *(-2.0*atan(rsph/DiskGravityDarkMatterR/Mpc_cm) + 
              2.0*log(1.0+rsph/DiskGravityDarkMatterR/Mpc_cm) +
              log(1.0+POW(rsph/DiskGravityDarkMatterR/Mpc_cm,2)));
  FtotR = -1*halo_g_of_r(rsph)/rsph * rcyl*LengthUnits;
  FtotR += -GravConst*DiskGravityStellarDiskMass*SolarMass*rcyl*LengthUnits
            /sqrt(POW(POW(rcyl*LengthUnits,2) + 
            POW(DiskGravityStellarDiskScaleHeightR*Mpc_cm +
                sqrt(POW(z,2) + 
                POW(DiskGravityStellarDiskScaleHeightz*Mpc_cm,2)),
                2),
            3));
  FtotR += -GravConst*DiskGravityStellarBulgeMass*SolarMass
           /POW(sqrt(POW(z,2) + POW(rcyl*LengthUnits,2)) + 
            DiskGravityStellarBulgeR*Mpc_cm,2)*rcyl*LengthUnits/
            sqrt(POW(z,2) +POW(rcyl*LengthUnits,2));

  /* Some error checking. */

  if (temperature < 0.0) 
    fprintf(stderr,"G_GSIG: temp = %"FSYM", P = %"FSYM", z = %"FSYM", zicm = %"FSYM", zicmf=%"FSYM", zsmall=%"FSYM", rcyl = %"FSYM"\n", 
      temperature, Pressure, z/Mpc_cm, zicm/Mpc_cm, zicmf, zsmall, rcyl*LengthUnits/Mpc_cm);
  if ((FtotR - FdPdR) > 0.0) { 
    fprintf(stderr,"G_GSIG: FtotR = %"FSYM", FdPdR = %"FSYM", P = %"FSYM",P2 = %"FSYM", Picm = %"FSYM", dr = %"FSYM", rcyl = %"FSYM", z = %"FSYM"\n", 
      FtotR, FdPdR, Pressure, Pressure2, Picm, r2-rcyl*LengthUnits, rcyl*LengthUnits/Mpc_cm, z/Mpc_cm);
    FdPdR = 0.0;
  } // end FtotR - FdPdr > 0.0 if

  /* Find circular velocity by balancing FG and dPdR against centrifugal force */

  vrot=sqrt(-rcyl*LengthUnits*(FtotR-FdPdR));

  if (denuse == densicm) vrot = 0.0;

  rotvel = vrot/VelocityUnits; //code units

} // end DiskForceBalance

// *************************************************************
// The functions integrated by qromb (parameter must be external)

// Given the in place radius and height, this returns density times
//    the stellar bulge force

double PbulgeComp_general(double rvalue, double zint)
{
  return (-MgasScale*SolarMass/
    (2*pi*POW(gScaleHeightR*Mpc_cm,2)*gScaleHeightz*Mpc_cm)*0.25/
    cosh(rvalue/gScaleHeightR/Mpc_cm) / cosh(fabs(zint)/gScaleHeightz/Mpc_cm)*
    GravConst*DiskGravityStellarBulgeMass*SolarMass/
    POW((sqrt(POW(zint,2) + POW(rvalue,2)) + 
         DiskGravityStellarBulgeR*Mpc_cm),2)*fabs(zint)/
    sqrt(POW(zint,2)+POW(rvalue,2)));
}

// Stellar Bulge functions

double PbulgeComp1(double zint)
{
  extern FLOAT rcyl;
  return PbulgeComp_general(rcyl*LengthUnits, zint);
}

double PbulgeComp2(double zint)
{
  extern FLOAT r2;
  return PbulgeComp_general(r2, zint);
}

// Given the in place radius and height, this returns density times
//    the stellar disk force

double PstellarComp_general(double rvalue, double zint)
{
  return (-MgasScale*SolarMass/
    (2*pi*POW(gScaleHeightR*Mpc_cm,2)*gScaleHeightz*Mpc_cm)*0.25/
    cosh(rvalue/gScaleHeightR/Mpc_cm) / cosh(fabs(zint)/gScaleHeightz/Mpc_cm)*
    GravConst*DiskGravityStellarDiskMass*SolarMass*
    (DiskGravityStellarDiskScaleHeightR*Mpc_cm + 
     sqrt(POW(zint,2) + POW(DiskGravityStellarDiskScaleHeightz*Mpc_cm,2)))*
    fabs(zint)/
    sqrt(POW(POW(rvalue,2) + 
       POW((DiskGravityStellarDiskScaleHeightR*Mpc_cm + 
      sqrt(POW(zint,2) + 
           POW(DiskGravityStellarDiskScaleHeightz*Mpc_cm,2)))
           ,2)
       ,3))/
    sqrt(POW(zint,2)+POW(DiskGravityStellarDiskScaleHeightz*Mpc_cm,2)));
}

// Stellar Disk functions

double PstellarComp1(double zint)
{
  extern FLOAT rcyl;
  return PstellarComp_general(rcyl*LengthUnits, zint);
}

double PstellarComp2(double zint)
{
  extern FLOAT r2;
  return PstellarComp_general(r2, zint);
}

double PDMComp_general(double rvalue, double zint){
/* --------------------------------------------------------
 * PDMComp_general
 * --------------------------------------------------------
 * General function for computing the DM contribution to
 * local (vertical) pressure on the gas in the galaxy's disk.
 * This returns the gas density at a given position
 * times the force on the gas due to the dark matter
 * --------------------------------------------------------- */

  float gas_density;
  float F;             // dark matter force
  float rsph; // 3D, spherical radius


  /* compute gas density */
  gas_density  = MgasScale*SolarMass / (8.0 * pi * POW(gScaleHeightR*Mpc_cm,2)*gScaleHeightz*Mpc_cm);
  gas_density /= (cosh(rvalue/gScaleHeightR/Mpc_cm)*cosh(fabs(zint)/gScaleHeightz/Mpc_cm));

  rsph = sqrt(rvalue*rvalue + zint*zint);

  /* fabs(zint) is because this is the force in the direction downward */
  F = -1*fabs(zint) * halo_g_of_r(rsph)/rsph;

  return gas_density * F;
}

/* DM pressure integration */
double PDMComp1(double zint){
  extern FLOAT rcyl;
  return PDMComp_general(rcyl*LengthUnits, zint);
}

double PDMComp2(double zint){
  extern FLOAT r2;
  return PDMComp_general(r2, zint);
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

/* Main integration routine called by DiskForceBalance to find Pressure */
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
  fprintf(stderr,"\t>> rcyl = %"FSYM", z = %"FSYM", z_icm = %"FSYM"\n", rcyl*LengthUnits/Mpc_cm, a/Mpc_cm, b/Mpc_cm);
  fprintf(stderr,"\t>> ss = %"FSYM", dss = %"FSYM"\n", ss, dss);
  ENZO_FAIL("FAILED IN QROMB IN GRID_GALAXYSIMULATIONINIALIZE\n");
  return -1.0;
}

/* -------------------- END of Functions for thermal pressure balance in the disk -------------------- */

/* --------------------Functions for initializing chemistry ------------------------------- */
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float& HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest)
{
  if (equilibrate) {
    /*  What temperature and density bins does the cell fall between? 
     *  'density' is in code units; 'temperature' is K
     *  Table arrays contain fractions.
     *  densities are returned in code units
     */
        
    // Start by assuming values are larger than those in table;
    // set to dim_size-1 for highest available value
    bool interpolate = true;
    int dens_indx, temp_indx, iter;
    dens_indx = temp_indx = EquilibriumTable.dim_size-1;
    
    for (iter=0; iter < EquilibriumTable.dim_size; ++iter) {
      if (density < EquilibriumTable.density[iter]) {
        dens_indx = iter-1;
        break;
      }
    }
    
    for (iter=0; iter<EquilibriumTable.dim_size; ++iter) {
      if (temperature < EquilibriumTable.temperature[iter]) {
        temp_indx = iter-1;
        break;
      }
    }

    // Density or temperature lower than in table
    if (dens_indx == -1) {
      dens_indx = 0;
      interpolate = false;
    }
    if (temp_indx == -1){
      temp_indx = 0;
      interpolate = false;
    }

    // Density or temp higher than table; unchanged from inital value
    if (dens_indx == EquilibriumTable.dim_size-1 ||
	      temp_indx == EquilibriumTable.dim_size-1)
      interpolate = false;

    if (interpolate) {
      HIdest = density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      
      DEdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      if (MultiSpecies > 1) {
	HMdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2Idest =  density*bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2IIdest =  density*bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
      if (MultiSpecies > 2) {
	DIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	DIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	HDIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
    } // end interpolate
    else { // don't interpolate; density and/or temp at edge of table
      HIdest =  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HIIdest = EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIdest = EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIdest = EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIIdest = EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      DEdest =  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx];

      if (MultiSpecies > 1) {
        HMdest =  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx];

        H2Idest = EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx];

        H2IIdest = EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
      if (MultiSpecies > 2) {
        DIdest =  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx];

        DIIdest = EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx];

        HDIdest = EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
    } // end no interpolation
  } // end if equilibrate
  else {
    HIdest = TestProblemData.HI_Fraction * density *
      TestProblemData.HydrogenFractionByMass;

    HeIdest = TestProblemData.HeI_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIdest = TestProblemData.HeII_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIIdest = (1.0 - TestProblemData.HydrogenFractionByMass) *
      density - HeIdest - HeIIdest;

    if(MultiSpecies > 1){
      HMdest = TestProblemData.HM_Fraction *
	            TestProblemData.HydrogenFractionByMass * density;
   
      H2Idest = 2 * TestProblemData.H2I_Fraction *
	            TestProblemData.HydrogenFractionByMass * density;
   
      H2IIdest = 2 * TestProblemData.H2II_Fraction 
	              * TestProblemData.HydrogenFractionByMass * density;
    }

    // HII density is calculated by subtracting off the various ionized fractions
    // from the total
    HIIdest = TestProblemData.HydrogenFractionByMass * density - HIdest;
    if (MultiSpecies > 1)
      HIIdest -= (HMdest + H2IIdest + H2Idest);

    // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
    // density for convenience) is calculated by summing up all of the ionized species.
    // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
    // calculating mass density, not number density (because the BaryonField values are 4x as
    // heavy for helium for a single electron)
    DEdest = HIIdest + 0.25*HeIIdest + 0.5*HeIIIdest;
    
    if (MultiSpecies > 1)
      DEdest += 0.5*H2IIdest - HMdest;
       
    DEdest = max(DEdest, tiny_number);
       
    // Set deuterium species (assumed to be a negligible fraction of the total, so not
    // counted in the conservation)
    if(MultiSpecies > 2){
      DIdest = HIdest * TestProblemData.DeuteriumToHydrogenRatio;
      DIIdest = HIIdest *	TestProblemData.DeuteriumToHydrogenRatio;
      HDIdest = H2Idest *	0.75 * TestProblemData.DeuteriumToHydrogenRatio;
    }

  } // end not equilibrate
}
double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2) {
    double interp;

    interp = f_x1y1*(x2-x)*(y2-y) + f_x2y1*(x-x1)*(y2-y) 
           + f_x1y2*(x2-x)*(y-y1) + f_x2y2*(x-x1)*(y-y1);

    interp *= 1/( (x2-x1)*(y2-y1) );

    return interp;

}

/* -------------------- BEGINNING OF Routines used for initializing the circumgalactic medium -------------------- */
/* 
   Computes halo gas density values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  Depending on the
   specific model chosen, different global parameters are needed (as set near the beginning
   of Grid::GalaxySimluationInitializeGrid).  Halo types are:

   GalaxySimulationGasHalo = 0  -- "zero CGM" - sets to a very low density/temperature
   GalaxySimulationGasHalo = 1  -- assuming hydrostatic equilibrium of CGM given an NFW dark matter halo 
                                   and a temperature as a function of radius set by the virial theorem.
   GalaxySimulationGasHalo = 2  -- assumes density, temperature set according to T = Tvir and entropy
                                   as a power-law function of radius.
   GalaxySimulationGasHalo = 3  -- as #2, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 4  -- assumes a hydrostatic equilibrium of CGM given an NFW dark matter halo
                                   and an entropy that is a power-law function of radius.
   GalaxySimulationGasHalo = 5  -- as #4, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 6  -- as #4, but the entropy distribution follows that for a precipitation-regulated NFW halo
                                   in Voit 2019 (ApJ)
   GalaxySimulationGasHalo = 7  -- previously reserved
   GalaxySimulationGasHalo = 8  -- a density and temperature profile fit to the entropy profiles in Voit 2019 (ApJ)     

   Inputs:  R - spherical radius, code units

   Returns:  density, grams/cm^3

   Note: using global variables w/following units:

   GalaxySimulationGasHalo: integer, >= 0
   GalaxySimulationGasHaloScaleRadius, units of Mpc_cm
   GalaxySimulationGasHaloDensity, units of grams/cm^3
   GalaxySimulationGasHaloTemperature, units of Kelvin
   GalaxySimulationGasHaloAlpha, power-law index; unitless
   GalaxySimulationGasHaloCoreEntropy, units of keV cm^2
   GalaxySimulationGasHaloMetallicity, units of Zsun
*/
double HaloGasDensity(FLOAT R, struct CGMdata& CGM_data){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low density */
   
    return densicm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets density assuming hydrostatic equilibrium using a temperature 
       as a function of radius given by virial theorem */
    
    double T0,haloDensity;
    T0 = HaloGasTemperature(GalaxySimulationGasHaloScaleRadius*Mpc_cm/LengthUnits,
			    CGM_data);
    haloDensity = GalaxySimulationGasHaloDensity*(T0/HaloGasTemperature(R, CGM_data));
    haloDensity /= POW((R*LengthUnits/GalaxySimulationGasHaloScaleRadius/Mpc_cm),3);
    return min(haloDensity,GalaxySimulationGasHaloDensity);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir, so
       n(r) = n_0 * (r/r_0)**(-alpha/(gamma-1))
       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */
    double scale_radius_cgs, this_radius_cgs, power_law_exponent;
    
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc_cm;
    this_radius_cgs = R*LengthUnits;
    power_law_exponent = -1.0*GalaxySimulationGasHaloAlpha/(Gamma-1.0);
    
    return GalaxySimulationGasHaloDensity*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent);
    
  } else if(GalaxySimulationGasHalo == 3){
    /* assumes entropy is a  power-law function of radius and T = Tvir that has a core (i.e., minimum 
       entropy value), so:

       n(r) = (Tvir / (Score + S_0*(r/r_0)^alpha))^(1/(gamma-1))

       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       Tvir is user-supplied temperature,
       Score is a user-supplied entropy,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */

    double scale_radius_cgs, this_radius_cgs, power_law_exponent, S_0, T_kev, n_0, this_number_density;

    // get radii in common set of units
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc_cm;
    this_radius_cgs = R*LengthUnits;
    power_law_exponent = -1.0*GalaxySimulationGasHaloAlpha/(Gamma-1.0);

    // now get number density using expression above

    T_kev = GalaxySimulationGasHaloTemperature*kboltzKeV;  // halo temperature in keV
    n_0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert n_0 to electron number density 
    S_0 = T_kev / POW(n_0,Gamma-1.0);   // S_0 in units of kev cm^2

    // get number density at this radius giving the requested info
    this_number_density = T_kev/(GalaxySimulationGasHaloCoreEntropy + S_0*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent));
    this_number_density = POW(this_number_density, 1.0/(Gamma-1.0));

    return this_number_density*mu*mh;  // return physical density
    
  } else if(GalaxySimulationGasHalo >= 4 && GalaxySimulationGasHalo <= 7){
    /* assumes entropy is a power-law function of radius OR a cored power-law function
       of radius and gas is in hydrostatic equilibrium w/the NFW halo.  */

    double this_radius_cgs, Rstart;
    int index;

    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr + 1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.n_rad[index]*mu*mh;  // return physical density

  } else if(GalaxySimulationGasHalo == 8){
    /* Eqn 24 in the Appendix of Voit 2019; a fit to the theoretical density profile of a precipitation-regulated NFW halo.
       Equation gives the electron number density, but pull the same trick as methods 2 & 3 and assume n_e = n */
    double this_radius_kpc, this_number_density;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    return this_number_density*mu*mh;  // return physical density
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasDensity().");
  }
  
} // end HaloGasDensity


/* 
   Computes halo gas temperature values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  The properties of each of the
   models are described immediately above this in the comments for the function HaloGasDensity(). 

   Inputs:  R - spherical radius, code units

   Returns:  Temperature, Kelvin
*/

double HaloGasTemperature(FLOAT R, struct CGMdata& CGM_data){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low temperature */
   
    return Ticm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets temperature as a function of radius given by virial theorem */

    return GravConst*DarkMatterMassEnclosed(R)*mu*mh/(3.0*kboltz*R*LengthUnits);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;
    
  } else if(GalaxySimulationGasHalo == 3){

    /* assumes entropy is a cored power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;

  } else if(GalaxySimulationGasHalo >= 4 && GalaxySimulationGasHalo <= 7){
    /* assumes entropy is a power-law function of radius and gas is in hydrostatic equilibrium */

    double this_radius_cgs;
    int index;
    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr+1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.T_rad[index];  // return temperature in Kelvin

  } else if(GalaxySimulationGasHalo == 8){
    /* Theoretical temperature profile of a precipitation-regulated NFW halo, using fits to n(r) and S(r) */
    double this_radius_kpc, this_number_density, this_entropy;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    this_entropy = GalaxySimulationGasHaloCoreEntropy * POW(this_radius_kpc, GalaxySimulationGasHaloAlpha);

    return this_entropy * POW(this_number_density, Gamma-1.0) / kboltzKeV; // units of K
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasTemperature().");
  }
  
}

/* Initializes arrays of number density, temperature, and radius for 
   choices of circumgalactic medium that require numerical integration based
   on user-defined parameters.  These quantities are stored in a global struct
   with arrays of size nbins for convenience (global within this file, at least).
   Rstop is the outer boundary of the integrtation in CGS; if negative, |Rstop|*R200 is used. 
   nbins defaults to 8192. */
 void halo_init(struct CGMdata& CGM_data, grid* Grid, FLOAT Rstop, int GasHalo_override){

  int halo_type=GalaxySimulationGasHalo;
  if (GasHalo_override) // not 0
    halo_type = GasHalo_override;
  
  if(halo_type < 4 || halo_type > 7) return;

  double k1, k2, k3, k4;
  double M, R200, rho_crit = 1.8788e-29*0.49;
  double Rstart;

  int index;
  
  M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
  
  R200 = pow(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  
  if (Rstop < 0)
    Rstop = fabs(Rstop)*R200;
  CGM_data.R_outer = Rstop;// integrate out to the virial radius of halo

  Rstart = 0.0; // you could force a different start if you liked

  // stepsize for RK4 integration and radial bins
  CGM_data.R_inner = Rstart;
  CGM_data.dr = (CGM_data.R_outer - CGM_data.R_inner)/ double(CGM_data.nbins); 
  
  if (halo_type < 6){
    double T0, n0, r0, dr;
    double this_n, this_radius, temperature;
    
    // set some quantities based on user inputs; this defines our integration
    T0 = GalaxySimulationGasHaloTemperature;
    n0 = GalaxySimulationGasHaloDensity / (mu*mh);
    r0 = GalaxySimulationGasHaloScaleRadius*Mpc_cm;
  
    // used for our numerical integration
    dr = CGM_data.dr;
    this_n = n0;
    this_radius = r0;

    // set the bin that we start at (otherwise it doesn't get set!)
    index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);
    CGM_data.n_rad[index] = this_n;
    CGM_data.T_rad[index] = T0;
    CGM_data.rad[index] = this_radius;

    /* starting at the point where the user has defined the radius, density, and 
       temperature, use RK4 to integrate the number density outward to R_outer using the expression 
       for dn_dr in another function.  Calculate the temperature using the entropy at this radius. */
    while(this_radius <= CGM_data.R_outer){
    
      // calculate RK4 coefficients.
      k1 = halo_dn_dr(this_radius,          this_n);
      k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1);
      k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2);
      k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3);

      // update density and radius
      this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
      this_radius += dr;  // new radius

      // calculate temperature at this radius using known entropy 
      temperature = halo_S_of_r(this_radius) * POW(this_n,Gamma-1.0);
      
      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);    
      CGM_data.n_rad[index] = this_n;
      CGM_data.T_rad[index] = temperature;
      CGM_data.rad[index] = this_radius;
    }
        
    /* now we do the same thing as above, but integrating inward to zero radius. */
    this_n = n0;
    this_radius = r0;
    dr *= -1.0;
    
    while(this_radius > CGM_data.R_inner){
      
      // calculate RK4 coefficients.
      k1 = halo_dn_dr(this_radius,          this_n);
      k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1);
      k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2);
      k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3);

      // update density and radius
      this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
      this_radius += dr;  // new radius
      
      // calculate temperature at this radius using known entropy 
      temperature = halo_S_of_r(this_radius) * POW(this_n,Gamma-1.0);

      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/(-1.0*dr) + 1.0e-3);

      if(index >= 0){
        CGM_data.n_rad[index] = this_n;
        CGM_data.T_rad[index] = temperature;
        CGM_data.rad[index] = this_radius;
      }
    }

  } else if (halo_type == 6) { 

    /* Integrate pressure assuming HSE & S(r) from Voit 2019, then convert to n & T,
       instead of integrating n(r) directly as with methods 4 & 5. This makes the boundary
       condition easier to handle.*/
    double dr, rmax, vcirc2_max;
    double this_press, this_ent, this_radius;//, this_n;
    double mu_ratio = 1.1/mu; // mu_e/mu
    double T_floor = 4e4; // IGM

    // boundary condition & quantities for integration
    dr = -1.0*CGM_data.dr;
    this_radius = R200;
    this_ent = halo_S_of_r(this_radius, Grid); // in erg*cm^2

    rmax = 2.163*R200/GalaxySimulationDMConcentration;
    vcirc2_max = GravConst * halo_mod_DMmass_at_r(rmax)/rmax;
    this_press = mu_ratio*POW(0.25*mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));

    // set the bin that we start at (otherwise it doesn't get set!)
    index = int((this_radius - CGM_data.R_inner)/(-1.0*dr) + 1.0e-3);
    CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma); // n_e ~ n_i
    CGM_data.T_rad[index] = POW( POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
    CGM_data.rad[index] = this_radius;

    // integrate inward from R200    
    while(this_radius > CGM_data.R_inner){
      
      // calculate RK4 coefficients.
      k1 = halo_dP_dr(this_radius,          this_press,             Grid);
      k2 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k1, Grid);
      k3 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k2, Grid);
      k4 = halo_dP_dr(this_radius + dr,     this_press + dr*k3,     Grid);
      
      // update radius, pressure, entropy
      this_radius += dr;  // new radius
      this_press += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4); // P @ new radius
      this_ent = halo_S_of_r(this_radius, Grid); // entropy @ new radius

      // store density and temperature in the struct
      index = int((this_radius - CGM_data.R_inner)/(-1.0*dr) + 1.0e-3);
      if (index >= 0){
        CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
        CGM_data.T_rad[index] = POW( POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
        CGM_data.rad[index] = this_radius;
      }
    }
    
    // Reset to boundary state
    dr = CGM_data.dr;
    this_radius = R200;
    this_ent = halo_S_of_r(this_radius, Grid); // in erg*cm^2

    rmax = 2.163*R200/GalaxySimulationDMConcentration;
    vcirc2_max = GravConst * halo_mod_DMmass_at_r(rmax)/rmax;
    this_press = mu_ratio*POW(0.25*mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));

    // Construct sigmoid to transition temperature to a constant
    double this_temp, this_dens;
    double deriv, r0, y0, y_offset, k;

    index = int((this_radius - CGM_data.R_inner)/(1.0*dr) + 1.0e-3);
    this_dens = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
    this_temp = POW( POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;

    deriv = (log10(this_temp) - log10(CGM_data.T_rad[index-1]))
          / (log10(this_radius) - log10(this_radius-dr));

    r0 = log10(this_radius);
    y0 = 2.0 * log10( T_floor / this_temp );
    assert (y0 < 0.0);
    y_offset = log10(this_temp) - y0/2.0;
    k = fabs(4.0/y0 * deriv);

    // Set constant dlog(P)/dlog(r)
    double prev_press, dlP_dlr, this_dPdr, press_vir;
    prev_press = mu_ratio * CGM_data.n_rad[index-1]/2.0 * kboltz*CGM_data.T_rad[index-1];
    press_vir = this_press;
    
    dlP_dlr = (log10(this_press) - log10(prev_press))
            / (log10(this_radius) - log10(this_radius-dr));
    assert (dlP_dlr < 0.0);
    
    while(this_radius <= CGM_data.R_outer){

      this_dPdr = this_press/this_radius * dlP_dlr;

      // update density and radius
      this_radius += dr;  // new radius
      this_dens = -2.0 * this_dPdr/(1.22*mh*halo_mod_g_of_r(this_radius)); // n_e = n_i
      this_temp = POW(10, sigmoid(log10(this_radius), r0, k, y0, y_offset));
      this_press = POW(10, dlP_dlr*log10(this_radius/R200) + log10(press_vir));

      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);    
      if (index < CGM_data.nbins) {
        CGM_data.n_rad[index] = this_dens;
        CGM_data.T_rad[index] = this_temp;
        CGM_data.rad[index] = this_radius;
      }
      else
        break;
    }
  }
    
  if (CGM_data.R_inner == 0) {
    // this integration acts a little squirrelly around r=0 because the mass values are garbage.  Cheap fix.
    CGM_data.rad[0]=CGM_data.rad[1];
    CGM_data.n_rad[0]=CGM_data.n_rad[1];
    CGM_data.T_rad[0]=CGM_data.T_rad[1];
  }
  
  return;
}

/* Halo entropy as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units.  output is entropy in CGS units (Kelvin cm^2) 
*/
double halo_S_of_r(FLOAT r){

  double Tvir, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc_cm;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4){

    return S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else if (GalaxySimulationGasHalo == 5){

    return Smin + S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else {
    ENZO_FAIL("halo_S_of_r: GalaxySimulationGasHalo set incorrectly.");
  }

}

/* More complex entropy profile from Voit 2019 that requires calculation of the cooling function.
   This one returns entropy in erg cm^2 instead of K cm^2 */
double halo_S_of_r(FLOAT r, grid* Grid){
  if (GalaxySimulationGasHalo == 6){

#ifdef USE_GRACKLE
    double M, C, r_vir, r_max, rho_crit = 1.8788e-29*0.49;
    double vcirc2, vcirc2_max;
    float Tgrav, Tgrav_therm;
    
    M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
    C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
    r_vir = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
    r_max = 2.163 * r_vir/C;
    
    vcirc2 = GravConst * halo_mod_DMmass_at_r(r) / r;
    vcirc2_max = GravConst * halo_mod_DMmass_at_r(r_max) / r_max;
    
    Tgrav = mu*mh * vcirc2 / kboltz; // 2x gravitational "temperature"
    Tgrav_therm = Tgrav / TemperatureUnits / ((Gamma-1.0)*mu); // code
  
    /* Calculate the cooling function Lambda using Grackle */
    float Lambda;
    float dens = mh/DensityUnits; // code
    float vx=0, vy=0, vz=0;
    float hi, hii, hei, heii, heiii, de, hm, h2i, h2ii, di, dii, hdi, metal; // species
    int dim=1;

    // setup_chem has densities in code, temperature in K
    setup_chem(dens, Tgrav, EquilibrateChem, de, hi, hii, hei, heii, heiii, hm, h2i, h2ii, di, dii, hdi);
    metal = GalaxySimulationGasHaloMetallicity * CoolData.SolarMetalFractionByMass * dens;

    // temporarily disable UV background; makes S(r) trend downward at large r instead of upward
    // because of low Tgrav
    int saved_UVB = grackle_data->UVbackground;
    grackle_data->UVbackground = 0;
    Grid->GrackleCustomCoolRate(1, &dim, &Lambda,
				&dens, &Tgrav_therm,
				&vx, &vy, &vz,
				&hi, &hii,
				&hei, &heii, &heiii,
				&de,
				&hm, &h2i, &h2ii,
				&di, &dii, &hdi,
				&metal);
    grackle_data->UVbackground = saved_UVB;

    // to cgs
    Lambda = fabs(Lambda) * POW(mh,2) * POW(LengthUnits,2) / ( POW(TimeUnits,3) * DensityUnits);

    /* Calculate entropy S(r) in erg cm^2 */
    double S_precip = POW(2*mu*mh, 1./3.) * POW(r*Lambda*GalaxySimulationGasHaloRatio/3.0, 2./3.);
    double S_nfw = 39. * vcirc2_max/1e10/4e4 * POW(r/r_vir, 1.1) / KEV_PER_ERG; // See Voit 19 Eqn 10 for assumptions

    // TODO blend with an entropy cap
    
    return (S_nfw + S_precip);
#else

    ENZO_FAIL("Grid_GalaxySimulationInitializeGrid: GasHalo 6 requires Grackle.");

#endif
  } else {
    ENZO_FAIL("halo_S_of_r: GalaxySimulationGasHalo set incorrectly.");
  }  
    
}

/* dEntropy/dr as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units; output is entropy gradient in CGS units (Kelvin*cm) */
double halo_dSdr(FLOAT r, double n){

  double Tvir, alpha, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc_cm;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4 || GalaxySimulationGasHalo == 5){

    // has units of Kelvin*cm, same deriv for both halo types (since constant drops out)
    return S0*GalaxySimulationGasHaloAlpha*
      POW(r/r0,GalaxySimulationGasHaloAlpha-1.0)/r0;

  } else {
    ENZO_FAIL("halo_dSdr: GalaxySimulationGasHalo set incorrectly.");
  }
}

/* a flexible, contortable sigmoid function. Gets used for log(K) */
double sigmoid(FLOAT x, FLOAT x0, double k, double y0, double y_off) {
  return y0 / (1 + exp(-k*(x-x0))) + y_off;
}

/* dn/dr as a function of radius and halo electron number density.  This quantity is calculated
   by assuming that gravity and pressure are in hydrostatic equilibrium in a halo with a specified 
   entropy profile S(r).

   Input is radius in CGM units and electron number density in units 
   of particles per cm^-3.  Output is dn/dr in CGS units, so particles per cm^4. 
*/
double halo_dn_dr(FLOAT r, double n){
  
  return -1.0*( n*1.22*mh*halo_g_of_r(r) + kboltz*POW(n,Gamma)*halo_dSdr(r,n) ) /
    ( Gamma * kboltz * halo_S_of_r(r) * POW(n, Gamma-1));
}

// double halo_dn_dr(double r, double n, double T) {
//   return -1.0 *  * 1.22*mh*n / (kboltz*T);
// }

double halo_dP_dr(FLOAT r, double P, grid* Grid) {
  return -1.0 * halo_mod_g_of_r(r) * 1.22*mh * POW( P/(1.1/mu) / halo_S_of_r(r,Grid),
						    1./Gamma );
}

/* halo gravitational acceleration as a function of radius.

   Input is the radius in CGS units and returns the MAGNITUDE of the 
   acceleration in CGS units.  */
double halo_g_of_r(FLOAT r){

  if (DiskGravityDarkMatterUseNFW)
    return GravConst*DarkMatterMassEnclosed(r)/(r*r);

  else if (DiskGravityDarkMatterUseB95) {
    double DMRad, DMDens;
    DMRad = DiskGravityDarkMatterR; // in Mpc
    DMDens = DiskGravityDarkMatterDensity; // in cgs

    return pi*GravConst*DMDens*POW(DMRad*Mpc_cm,3)/POW(r,2) // NOT r^3
            *(-2.0*atan(r/DMRad/Mpc_cm)
              +2.0*log(1.0+r/DMRad/Mpc_cm)
              +log(1.0+POW(r/DMRad/Mpc_cm,2))
              );
  }
}

double halo_mod_g_of_r(FLOAT r){
  return GravConst*halo_mod_DMmass_at_r(r)/(r*r);
}

/* A modified NFW halo, where small radii (r <= 2.163*Rs) have enclosed mass
   that satisfies v_circ(r) = v_circ,max = v_circ(2.163*Rs). Beyond this
   boundary radius, the enclosed mass is that of a standard NFW halo.
*/
double halo_mod_DMmass_at_r(FLOAT r){

  double M, C, R200, Rs, Rmax;
  double rho_crit = 1.8788e-29*0.49;
  
  M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
  C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
  
  R200 = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  Rs = R200/C;  // scale radius of NFW halo in CGS
  Rmax = 2.163*Rs;
  
  if (r <= Rmax) {
    return r / Rmax * DarkMatterMassEnclosed(Rmax);
  }
  else {
    return DarkMatterMassEnclosed(r);
  }
}

/* -------------------- END OF Routines used for initializing the circumgalactic medium -------------------- */


