/***********************************************************************
/
/  GRID CLASS (INITIALIZE ISOLATED DISK GALAXY)
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1:
/
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
float gasdev();

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

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.672e-8,
               pi = 3.14159, mh = 1.6726e-24, kboltz = 1.3807e-16;
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
    NFWMass[i] = 4.0*pi*HaloDensity[sphere]*DensityUnits*
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
    fprintf(fptr, "%d %"GOUTSYM" %g %g %g %g %g %g %g\n", i, NFWRadius[i]*LengthUnits, 
	 NFWDensity[i], NFWMass[i], NFWPressure[i], NFWTemp[i], NFWSigma[i],
         mean_overdensity, NFWVelc[i]);
    if (mean_overdensity > 199 && m200 == 0) {
      m200 = NFWMass[i]/SolarMass;
      rhomean = mean_overdensity;
      r_vir = NFWRadius[i]*LengthUnits;
      i_vir = i;
      if (i_vir != 0) {
	printf("Halo parameters wrong: mean_overdensity=%g\n",
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
  fprintf(fptr, "#m200 = %g, r_vir = %g, T_vir=%g, r_sphere=%g, vfac=%g,rhomean=%g\n", 
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
    
    printf("rho_sphere = %g, rho_p = %g\n", HaloDensity[0],
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
      printf("npart=%d\n", NumberOfParticles);
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
	      
	      FLOAT R = sqrt(xpos*xpos+ypos*ypos);

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
    printf("GalaxyDiskInitialize: NumberOfParticles = %d\n", 
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
