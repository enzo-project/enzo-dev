/***********************************************************************
/
/  GRID CLASS (INITIALIZE ISOLATED DISK GALAXY)
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1:  May 2010, Tom Abel added exponential disk within NFW halo
/
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
float gasdev();

double BESSI0(double X);
double BESSI1(double X);
double BESSK0(double X);
double BESSK1(double X);

static int CollapseTestParticleCount = 0;


int grid::GalaxyDiskInitializeGrid(int NumberOfHalos,
				   FLOAT HaloRadius[MAX_SPHERES],
				   FLOAT HaloCoreRadius[MAX_SPHERES],
				   float HaloDensity[MAX_SPHERES],
				   float HaloTemperature[MAX_SPHERES],
				   FLOAT HaloPosition[MAX_SPHERES][MAX_DIMENSION],
				   float HaloSpin[MAX_SPHERES],
				   float HaloVelocity[MAX_SPHERES][MAX_DIMENSION],
				   float HaloAngVel[MAX_SPHERES],
				   float HaloMagneticField,
				   FLOAT DiskRadius[MAX_SPHERES],
				   FLOAT DiskHeight[MAX_SPHERES],
				   float DiskDensity[MAX_SPHERES],
				   float DiskTemperature[MAX_SPHERES],
				   float DiskMassFraction[MAX_SPHERES],
				   float DiskFlaringParameter[MAX_SPHERES],
				   int   GalaxyType[MAX_SPHERES],
				   int   UseParticles, int UseGas,
				   float UniformVelocity[MAX_DIMENSION],
				   float MediumTemperature, float MediumDensity, int level)
{
  /* declarations */

  int dim, i, j, k, m, field, sphere, size;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* create fields */

  NumberOfBaryonFields = 0;
  int ivel;
  int phip_num;
  if (UseGas) {
    FieldType[NumberOfBaryonFields++] = Density;
    ivel = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    FieldType[NumberOfBaryonFields++] = Velocity2;
    FieldType[NumberOfBaryonFields++] = Velocity3;
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    if (HydroMethod == MHD_RK) {
      FieldType[NumberOfBaryonFields++] = Bfield1;
      FieldType[NumberOfBaryonFields++] = Bfield2;
      FieldType[NumberOfBaryonFields++] = Bfield3;
      FieldType[NumberOfBaryonFields++] = PhiField;
    }
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
  }

  if(UseDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
    FieldType[NumberOfBaryonFields++] = DebugField;  
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (UseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  /* Set various units. */

  double  mu = Mu; // assume fully ionized cosmic gas
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits*pow(LengthUnits,3);

  /* Compute NFW profile-info. The parameters are HaloCoreRadius (which is
     the "knee radius") and HaloDensity (the overdensity ar the knee
     radius).  The pressure is computed by integrated the hydro-static
     equation and from this comes the temperature and the dm velocity
     dispersion. */

#define NFW_POINTS 500
  double NFWMass[NFW_POINTS], NFWPressure[NFW_POINTS], NFWTemp[NFW_POINTS], x1,
        NFWDensity[NFW_POINTS], NFWSigma[NFW_POINTS], m200, NFWVelc[NFW_POINTS];
  FLOAT NFWRadius[NFW_POINTS];

  // Compute the parameter for circular velocity
  // Reference: Springel & White (1999)

  double l_spin = HaloSpin[0];
  double c = HaloRadius[0]/HaloCoreRadius[0];

  double g_c = 0, dx = c/1000.0;
  double xi;
  for (int i = 0; i < 1000; i++) {
    xi = (i+0.5)*dx;
    g_c += sqrt(log(1.0+xi) - xi/(1.0+xi))*pow(xi,1.5)/pow(1.0+xi,2);
  }
  g_c *= dx;
  double f_c = c*(0.5-0.5/pow(1.0+c,2)-log(1.0+c)/(1.0+c))/pow(log(1.0+c)-c/(1.0+c),2);
  double f_s = 1.5*l_spin*sqrt(2.0*c/f_c)/g_c*pow(log(1.0+c)-c/(1.0+c),1.5);
  double delta_c = 200.0/3.0*pow(c,3)/(log(1.0+c)-c/(1.0+c));

  double dpdr = 0, dpdr_old;
  double r_vir = 0;
  double mean_overdensity;
  int i_vir;
  sphere = 0;
  m200   = 0;
  //NFWPressure[0] = MediumDensity * DensityUnits * kboltz * MediumTemperature / (mu * mh);
  FILE *fptr = fopen("NFWProfile.out", "w");
  float rhomean;
  for (i = 0; i < NFW_POINTS; i++) {
    NFWRadius[i] = HaloRadius[sphere]*pow(10, -3*(float(i)/NFW_POINTS));
    x1 = NFWRadius[i]/HaloCoreRadius[sphere];
    NFWDensity[i] = HaloDensity[sphere]*DensityUnits/(x1*(1.0+x1)*(1.0+x1));
    NFWMass[i] = 4.0*M_PI*HaloDensity[sphere]*DensityUnits*
		pow(HaloCoreRadius[sphere]*LengthUnits, 3) *
		(log(1.0+x1) - x1/(x1+1.0));  // in g
    dpdr_old = dpdr;
    dpdr = GravConst * NFWMass[i] * NFWDensity[i] /
           pow(NFWRadius[i]*LengthUnits, 2);
    if (i == 0) {
      NFWPressure[0] = NFWDensity[0]*kboltz*MediumTemperature/(mu*mh);
      //NFWPressure[0] = 0.5*NFWDensity[0]*GravConst*NFWMass[0]/(NFWRadius[0]*LengthUnits);
    }
    if (i > 0) {
      NFWPressure[i] = NFWPressure[i-1] -
	0.5*(dpdr+dpdr_old)*(NFWRadius[i]-NFWRadius[i-1])*LengthUnits;
    }
    NFWTemp[i] = NFWPressure[i]*mu*mh/(kboltz*NFWDensity[i]); // K
    NFWSigma[i] = sqrt(NFWPressure[i]/NFWDensity[i]); // cm/s
      //sqrt(kboltz * NFWTemp[i] / (mu * mh));  // in cm/s

    NFWVelc[i] = f_s*sqrt(GravConst*NFWMass[i]/(NFWRadius[i]*LengthUnits));
    
    mean_overdensity = 3.0*delta_c/pow(x1,3)*(log(1.0+x1)-x1/(x1+1.0));
    fprintf(fptr, "%"ISYM" %"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", i, NFWRadius[i]*LengthUnits, 
	 NFWDensity[i], NFWMass[i], NFWPressure[i], NFWTemp[i], NFWSigma[i],
         mean_overdensity, NFWVelc[i]);
    if (mean_overdensity > 199 && m200 == 0) {
      m200 = NFWMass[i]/SolarMass;
      rhomean = mean_overdensity;
      r_vir = NFWRadius[i]*LengthUnits;
      i_vir = i;
      if (i_vir != 0) {
	printf("Halo parameters wrong: mean_overdensity=%"GSYM"\n",
	       mean_overdensity);
	return FAIL;
      }
    }
  }

 

  /* Find the average density */
  double rho_nfw = 3.0*NFWMass[i_vir]/(4.0*M_PI*pow(NFWRadius[i_vir]*LengthUnits,3));

  /* Renormalize velocity using virial relation */
  double Ek = 0.0;
  for (int i = 1; i < NFW_POINTS; i++) {
    Ek += 2.0*M_PI*NFWDensity[i]*pow(NFWSigma[i],2)*pow(NFWRadius[i],2)*
      (NFWRadius[i-1]-NFWRadius[i])*pow(LengthUnits,3);
  }
  
  double vfac = sqrt(0.5*GravConst*NFWMass[0]*NFWMass[0]*f_c/r_vir/Ek);
  //double vfac = sqrt(GravConst*m200*SolarMass/r_vir)/NFWSigma[0];
  for (int i = 0; i < NFW_POINTS; i++) {
    NFWSigma[i] *= vfac;
  }
  double T_vir = 0.5*mh/kboltz*GravConst*m200*SolarMass/r_vir; 
  fprintf(fptr, "#m200 = %"GSYM", r_vir = %"GSYM", T_vir=%"GSYM", r_sphere=%"GSYM", vfac=%"GSYM",rhomean=%"GSYM"\n", 
	  m200, r_vir, T_vir, HaloRadius[0]*LengthUnits, vfac, rhomean);
  fclose(fptr);

  float CosmologySimulationInitialFractionHII   = 1.2e-5;
  float CosmologySimulationInitialFractionHeII  = 1.0e-14;
  float CosmologySimulationInitialFractionHeIII = 1.0e-17;
  float CosmologySimulationInitialFractionHM    = 2.0e-9;
  float CosmologySimulationInitialFractionH2I   = 2.0e-20;
  float CosmologySimulationInitialFractionH2II  = 3.0e-14;

  /* Loop over the set-up twice, once to count the particles, the second
     time to initialize them. */

  int SetupLoopCount, npart = 0;
  int UseBH = 0; //= (level == 0) ? 1 : 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(UseParticles+UseBH, 1);
       SetupLoopCount++) {

    /* Set densities */
    
    double ParticleCount = 0;
    FLOAT dx = CellWidth[0][0]*LengthUnits;
    double ParticleMeanDensity = NFWMass[i_vir]/1e4/pow(dx,3);
    ParticleMeanDensity /= DensityUnits;
    
    printf("rho_sphere = %"GSYM", rho_p = %"GSYM"\n", HaloDensity[0],
	   ParticleMeanDensity);
    
    /* Set particles. */
    
    if ((UseParticles > 0 || UseBH) && SetupLoopCount > 0) {
      
      /* If particles already exist (coarse particles), then delete. */
      
      if (NumberOfParticles > 0)
	this->DeleteParticles();

      /* Use count from previous loop to set particle number. */

      NumberOfParticles = npart+UseBH;
      npart = 0;
      
      /* Allocate space. */
      printf("npart=%"ISYM"\n", NumberOfParticles);
      this->AllocateNewParticles(NumberOfParticles);
      
      /* Particle values will be set below. */
      
    } // end: particle initialization

    /* Set up the baryon field. */
    
    size = 1;
    for (dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
    }
    
    if (SetupLoopCount == 0 && UseGas) {
      for (field = 0; field < NumberOfBaryonFields; field++) {
	if (BaryonField[field] == NULL) {
	  BaryonField[field] = new float[size];
	}
      }
    }
    /* Loop over the mesh. */
    
    float density, Velocity[MAX_DIMENSION], temperature, sigma;
    FLOAT r, x, y = 0, z = 0;
    int n = 0;
    float f_b = 1.0/10.0;
    float RotVelocity[3];
    FLOAT xpos, ypos, zpos, drad, cosphi, sinphi, sintheta; 

    
    // this sets up an isothermal disk in equilbrium vertically and radially within a static NFW halo
    // Central gas surface density Sigma_0 using DiskRadius for DiskScale Radius
    // from M_disk = 6 Pi Sigma_0 R_D^2  from  Sigma(R) = 2 rho_0 h_0 (1 + R/R_D/F) Exp[ -R/R_D ]
    // where F is the DiskFlaringParameter; 
    // rho_0 is central midplane density; h_0 is central scale height and R_D is the disk scale radius
    // R: is the cylindrical radius not spherical
    
    double R_D = DiskRadius[sphere] * LengthUnits;
    double R;
    double F = DiskFlaringParameter[sphere];
    double Sigma_0 = (double)DiskMassFraction[sphere]*NFWMass[0]/(R_D*R_D*2*M_PI)*F/(2+F);
    double c_s;
    double tgamma = Gamma;
    if (EOSType == 3)
      tgamma = 1.;

    c_s = sqrt(tgamma/Mu/mh * kboltz * DiskTemperature[0]); // simple gamma law (gamma should be close to one to be in equilibrium)
    double rho_0 = Sigma_0*Sigma_0*M_PI*GravConst/(c_s*c_s);
    double h_0   = Sigma_0/2/rho_0;

    double rhoc, nd, SigmaR, Omega, Mdr, hr, Vrot, ToomreQ, yd, vphidisk, vdm, vpress;
    if (GalaxyType[sphere] == 5) { 
      fptr = fopen("DiskProfile.out", "w");
      printf("#mdisk= %g, r_disk = %g kpc, T_disk=%g K, central height h_0=%g kpc, \nrho_0=%g [g/cm3], c_s=%g km/s, Sigma_0=%g Msun/pc^2\n", 
	     DiskMassFraction[sphere]*NFWMass[0]/SolarMass, 
	     R_D/kpc, DiskTemperature[0], h_0/kpc, rho_0, c_s/1e5, Sigma_0/SolarMass*pc*pc); 
      fprintf(fptr, "#mdisk= %g, r_disk = %g kpc, T_disk=%g K, central height h_0=%g kpc, \nrho_0=%g [g/cm3], c_s=%g km/s, Sigma_0=%g Msun/pc^2\n", 
	     DiskMassFraction[sphere]*NFWMass[0]/SolarMass, 
	     R_D/kpc, DiskTemperature[0], h_0/kpc, rho_0, c_s/1e5, Sigma_0/SolarMass*pc*pc); 

      fprintf(fptr, "   Radius [kpc]  Sigma(R) [Msun/pc^2]  n(R,z=0) [cm^{-3}]  M_disk(<R) [Msun] M_tot(<R) [Msun]  h(r) [pc] V_rot [km/s] V_dm [km/s] V_disk [km/s]  Vpress [km/s]  Toomre-Q    \n" );

      for (i = 0; i < NFW_POINTS; i++) {
	R = HaloRadius[sphere]*pow(10, -3*(float(i)/NFW_POINTS))*LengthUnits;
	SigmaR = Sigma_0*exp(-R/R_D)*(1+R/R_D/F);
	rhoc = 2*M_PI*GravConst*SigmaR*SigmaR/c_s/c_s;
	hr = SigmaR/2/rhoc;
	Mdr = M_PI*2*Sigma_0/F *( (exp(R/R_D)-1)*(2+F)*R_D*R_D - (2+F)*R*R_D - R*R)*exp(-R/R_D);
	yd = R/R_D;
	vphidisk = sqrt(4.*M_PI*GravConst*Sigma_0*R_D*yd*yd
			*(BESSI0(yd)*BESSK0(yd)-BESSI1(yd)*BESSK1(yd)));
	vpress = sqrt(2*R*(R+(F-1)*R_D)/(R_D*(R+F*R_D)))*c_s;
	vdm =  sqrt(GravConst*NFWMass[i]/R) ;
	Vrot = max(vdm + vphidisk - vpress,0); // do not go below 0
	Omega = Vrot/2/M_PI/R;
	ToomreQ = c_s*Omega/SigmaR/GravConst;
	fprintf(fptr, "%"ISYM" %"GOUTSYM"\t %g   \t %g\t  %g\t  %g\t  %g\t       %g\t  %g\t   %g\t %g \t %g\n", 
		i, R/kpc, SigmaR/SolarMass*pc*pc, rhoc/Mu/mh, Mdr/SolarMass, (NFWMass[i]+Mdr)/SolarMass, 
		hr/pc, Vrot/1e5, vdm/1e5,  vphidisk/1e5, -vpress/1e5,ToomreQ
		  );
	
      }	
      fclose(fptr);

    }


    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	for (i = 0; i < GridDimension[0]; i++, n++) {
	  
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	  density = MediumDensity;
	  temperature = MediumTemperature;
	  sigma = 0;
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    Velocity[dim] = 0;
	  }

	  for (sphere = 0; sphere < NumberOfHalos; sphere++) {
	    
	    /* Find distance from center. */
	    
	    r = sqrt(pow(fabs(x-HaloPosition[sphere][0]), 2) +
		     pow(fabs(y-HaloPosition[sphere][1]), 2) +
		     pow(fabs(z-HaloPosition[sphere][2]), 2) );
	    r = max(r, 0.1*CellWidth[0][0]);
	    
	    if (r < HaloRadius[sphere]) {
	      xpos = x-HaloPosition[sphere][0];
	      ypos = y-HaloPosition[sphere][1];
	      zpos = z-HaloPosition[sphere][2];
	      
	      R = sqrt(xpos*xpos+ypos*ypos);
	      //	      R = max(R, 0.1*CellWidth[0][0]);

	      // compute the azimuthal angle
	      cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
	      sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);

	      // compute the polar angle
	      sintheta = sqrt(xpos*xpos+ypos*ypos)/sqrt(xpos*xpos+ypos*ypos+zpos*zpos);

	      /* 1) Uniform */
	      
	      if (GalaxyType[sphere] == 1) {
		density = HaloDensity[sphere];
		temperature = HaloTemperature[sphere];
	      }

	      /* 2) r^-2 power law */
	      
	      if (GalaxyType[sphere] == 2) {
		density = HaloDensity[sphere]*pow(r/HaloRadius[sphere], -2);
		temperature = HaloTemperature[sphere];
	      }
	      
	      /* 3) NFW profile (use look-up table for temperature and
		 velocity dispersion)*/
	      
	      if (GalaxyType[sphere] == 3) {
		double vphi, vphim, vphip;
		x1 = r/HaloCoreRadius[sphere];
		density = f_b*HaloDensity[sphere]/((x1+0.2)*(1.0+x1)*(1.0+x1));
		for (m = 1; m < NFW_POINTS; m++) {
		  if (r >= NFWRadius[m]) {
		    /* Interpolate temperature, velocity dispersion
		       circular velocity */
		    //temperature = NFWTemp[m] + (NFWTemp[m-1] - NFWTemp[m])*
		    //(r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		    //temperature = 2e4;
		    temperature = T_vir;
		    sigma = NFWSigma[m] + (NFWSigma[m-1] - NFWSigma[m])*
		      (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		    vphi = NFWVelc[m] + (NFWVelc[m-1] - NFWVelc[m])*
		      (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		    vphi /= VelocityUnits;
		    break;
		  }
		}
		RotVelocity[0] = -vphi*sinphi*sintheta;
		RotVelocity[1] = vphi*cosphi*sintheta;
		Velocity[0] = RotVelocity[0];
		Velocity[1] = RotVelocity[1];
		// Add random velocity dispersion
		if (sigma != 0) {
		  for (dim = 0; dim < GridRank; dim++) {
		    Velocity[dim] +=
		      0.5*gasdev()*sigma/VelocityUnits;
		  }
		}
		
	      }	    

	      /* 4) Uniform disk embeded in NFW profile */
	      if (GalaxyType[sphere] == 4 && 
		  R < DiskRadius[sphere] &&
		  fabs(zpos) < DiskHeight[sphere]) {
		double vphi, vphim, vphip, m_disk;
		density = DiskDensity[sphere];		
		temperature = DiskTemperature[sphere];
		// add random perturbation to density and temperature
		density *= (1.0+0.01*gasdev());
		temperature *= (1.0+0.01*gasdev());
		for (m = 1; m < NFW_POINTS; m++) {
		  if (R >= NFWRadius[m]) {
		    m_disk = DiskDensity[sphere]*M_PI*R*R*2.0*DiskHeight[sphere];
		    m_disk *= MassUnits;
		    vphi = sqrt(GravConst*(NFWMass[m]+m_disk)/(R*LengthUnits));
		    vphi /= VelocityUnits;
		    break;
		  }
		}
		Velocity[0] = -vphi*sinphi;
		Velocity[1] = vphi*cosphi;
	      }	    

	      /* 5) Exponential disk in NFW halo */

	      if (GalaxyType[sphere] == 5 && R < HaloRadius[sphere]) { 
		// this sets up an isothermal disk in equilbrium vertically and radially within a static NFW halo
		// Central gas surface density Sigma_0 using DiskRadius for DiskScale Radius
		// from M_disk = 6 Pi Sigma_0 R_D^2  from  Sigma(R) = 2 rho_0 h_0 (1 + R/R_D) Exp[ -R/R_D ]
		// rho_0 is central midplane density; h_0 is central scale height and R_D is the disk scale radius
		// R: is the cylindrical radius not spherical

		R_D = DiskRadius[sphere];
		SigmaR = Sigma_0*exp(-R/R_D)*(1+R/R_D/F);
		rhoc = 2*M_PI*GravConst*SigmaR*SigmaR/c_s/c_s;
		hr = SigmaR/2/rhoc;

		// if the central density in the midplane is not linearly decreasing with radius
		// you do not have enough resolution for your parameter choice
		density += rhoc
		  /pow(cosh(zpos/(hr/LengthUnits)), 2)/DensityUnits;
		  //		  /pow(cosh(zpos/(hr/LengthUnits)), 2)/DensityUnits;

		if ((r < CellWidth[0][0]) && (CellWidth[0][0] > R_D)) 
		    density = 2*MinimumMassForRefinement[0]/CellWidth[0][0]/CellWidth[0][0]/CellWidth[0][0];


		// rotational velocity added by disk:
		yd = R/DiskRadius[0];
		vphidisk = sqrt(4.*M_PI*GravConst*Sigma_0*R_D*LengthUnits*yd*yd
				       *(BESSI0(yd)*BESSK0(yd)-BESSI1(yd)*BESSK1(yd)));
		// printf("%g : %g %g %g %g \n", yd,BESSI0(yd),BESSK0(yd),BESSI1(yd),BESSK1(yd)); // works!

		
		// pressure gradient subtracts some from the rotational velocity
		    vpress = sqrt(2*R*(R+(F-1)*R_D)/(R_D*(R+F*R_D)))*c_s;

		double vphi=0.;
		temperature = DiskTemperature[0];
		// add random perturbation to density and temperature
		//		density *= (1.0+0.01*gasdev());
		//		temperature *= (1.0+0.01*gasdev());
		for (m = 1; m < NFW_POINTS; m++) {
		  if (R >= NFWRadius[m]) {
		    vphi = sqrt(GravConst*(NFWMass[m])/(R*LengthUnits));
		    // printf("R: %g [kpc] : %g %g %g \n", R*LengthUnits/kpc, vphi, vphidisk, vphiP);
		    break;
		  }
		}

		    vphi += vphidisk; //  disk mass 
		    //		    vphi -= vpress;    // - pressure gradient
		
		vphi /= VelocityUnits;

		Velocity[0] = -vphi*sinphi;
		Velocity[1] = vphi*cosphi;

	      }

	    } // end: if (r < HaloRadius)
	  } // end: loop over spheres

	  /* Set BaryonField if using gas */

	  if (UseGas) {
	    
	    /* Set density. */
	    
	    BaryonField[iden][n] = density;
	    
	    /* If doing multi-species (HI, etc.), set these. */
	    
	    if (MultiSpecies > 0) {
	      
	      BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
		CoolData.HydrogenFractionByMass * BaryonField[0][n];
	      BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
		BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	      BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
		BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	      BaryonField[HeINum][n] = 
		(1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
		BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
	      
	      if (MultiSpecies > 1) {
		BaryonField[HMNum][n] = CosmologySimulationInitialFractionHM*
		  BaryonField[HIINum][n]* pow(temperature,float(0.88));
		BaryonField[H2IINum][n] = CosmologySimulationInitialFractionH2II*
		  2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
		BaryonField[H2INum][n] = CosmologySimulationInitialFractionH2I*
		  BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1);
	      }
	      
	      BaryonField[HINum][n] = 
		CoolData.HydrogenFractionByMass*BaryonField[0][n]
		- BaryonField[HIINum][n];
	      if (MultiSpecies > 1)
		BaryonField[HINum][n] -= BaryonField[HMNum][n]
		  + BaryonField[H2IINum][n]
		  + BaryonField[H2INum][n];
	      
	      BaryonField[DeNum][n] = BaryonField[HIINum][n] + 
		0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	      if (MultiSpecies > 1)
		BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - 
		  BaryonField[HMNum][n];
	      
	      /* Set Deuterium species (assumed to be negligible). */
	      
	      if (MultiSpecies > 2) {
		BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
		  BaryonField[HINum][n];
		BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
		  BaryonField[HIINum][n];
		BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
		  BaryonField[H2INum][n];
	      }
	    }

	    /* Set Velocities. */
	  
	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[ivel+dim][n] = Velocity[dim] + UniformVelocity[dim];
	    
	    /* Set energy (thermal and then total if necessary). */
	  
	    BaryonField[ietot][n] = temperature/TemperatureUnits/
	      ((Gamma-1.0)*mu);

	    if (DualEnergyFormalism)
	      BaryonField[ieint][n] = BaryonField[ietot][n];

	    for (dim = 0; dim < GridRank; dim++)
	      BaryonField[ietot][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);

	    if (HydroMethod == MHD_RK) {
	      BaryonField[iBx][n] = 0.0;
	      BaryonField[iBy][n] = 0.0; //HaloMagneticField;
	      BaryonField[iBz][n] = HaloMagneticField;
	      BaryonField[iPhi][n] = 0.0;
	      BaryonField[ietot][n] += 0.5*HaloMagneticField*HaloMagneticField/density;
	    }
	  } // if (UseGas)
	  
	  /* Set particles if being used (generate a number of particle
	     proportional to density). */

	  if (UseParticles) {
	    for (sphere = 0; sphere < NumberOfHalos; sphere++) {
	      
	      /* Find distance from center. */
	      
	      r = sqrt(pow(fabs(x-HaloPosition[sphere][0]), 2) +
		       pow(fabs(y-HaloPosition[sphere][1]), 2) +
		       pow(fabs(z-HaloPosition[sphere][2]), 2) );
	      r = max(r, 0.1*CellWidth[0][0]);

	      if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
		  j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
		  k >= GridStartIndex[2] && k <= GridEndIndex[2] &&
		  r < HaloRadius[sphere]) {
		
		//ParticleCount += density/pow(float(RefineBy), GridRank*level)/f_b;
		ParticleCount += density/f_b/ParticleMeanDensity;
		while (ParticleCount > 1) {
		  if (SetupLoopCount > 0) {
		    ParticleMass[npart] = ParticleMeanDensity;
		    //*pow(float(RefineBy), GridRank*level);
		    ParticleNumber[npart] = CollapseTestParticleCount++;
		    ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
		    
		    /* Set random position within cell. */
		    
		    ParticlePosition[0][npart] = x + 
		      CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		    ParticlePosition[1][npart] = y +
		      CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		    ParticlePosition[2][npart] = z +
		      CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		    
		    /* Set bulk velocity. */
		    
		    for (dim = 0; dim < GridRank; dim++)
		      ParticleVelocity[dim][npart] = 
			RotVelocity[dim]+UniformVelocity[dim];
		    
		    /* Add random velocity; */
		    
		    if (sigma != 0) {
		      for (dim = 0; dim < GridRank; dim++) {
			ParticleVelocity[dim][npart] += 
			  gasdev()*sigma/VelocityUnits;
		      }
		    }
		  } // if (SetupLoopCount>0)
		  npart++;
		  ParticleCount -= 1.0;
		} // end: while

	      } // end: if (i >= GridStarIndex...
	    } // for sphere
	  } // end: if (UseParticles...
	} // end for i
      } // end for j
    } // end for k
  } // end loop SetupLoopCount

  if (UseParticles) {
    printf("GalaxyDiskInitialize: NumberOfParticles = %"ISYM"\n", 
	   NumberOfParticles);
  }
  
  if (level == 0 && UseBH) {
    double m = 1e8;
    m *= 1.989e33;
    m /= MassUnits;
    ParticleType[NumberOfParticles-1] = PARTICLE_TYPE_DARK_MATTER;
    ParticleNumber[NumberOfParticles-1] = CollapseTestParticleCount++;
    ParticleMass[NumberOfParticles-1] = m/pow(CellWidth[0][0],3);
    ParticlePosition[0][NumberOfParticles-1] = 0.5;
    ParticlePosition[1][NumberOfParticles-1] = 0.5;
    ParticlePosition[2][NumberOfParticles-1] = 0.5;
    ParticleVelocity[0][NumberOfParticles-1] = 1e6/VelocityUnits;
    ParticleVelocity[1][NumberOfParticles-1] = 0.0;
    ParticleVelocity[2][NumberOfParticles-1] = 0.0;
  }
  
  return SUCCESS;
}
