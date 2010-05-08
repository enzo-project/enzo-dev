/***********************************************************************
/
/  GRID CLASS (ADD THE CONTENTS OF THIS GRID TO THE GIVEN RADIAL PROFILE)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/Fluxes.h"
#include "../enzo/GridList.h"
#include "../enzo/ExternalBoundary.h"
#include "../enzo/Grid.h"
#include "../enzo/CosmologyParameters.h"
#include "../enzo/units.h"

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int field, int farray[], int numfields);


int grid::AddToRadialProfile(FLOAT SphereCenter[MAX_DIMENSION], 
			     float SphereRadius,
			     FLOAT MeanVelocity[MAX_DIMENSION][3],
			     int NumberOfBins,
			     FLOAT ProfileRadius[], 
			     FLOAT ProfileValue[][MAX_PROFILES],
			     FLOAT ProfileWeight[][MAX_PROFILES],
			     char *ProfileName[MAX_PROFILES],
			     AnalyzeClusterParameters *parameters)
{

  /* Return if the data is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int i, j, k, index, m, n, dim;
  int ip, im, jp, jm, kp, km;
  FLOAT radius2, delx, dely, delz, InertiaWeight, radial_vel, AngWeight,
        velx, vely, velz, vel[MAX_DIMENSION], mu, number_density;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, kboltz = 1.38e-16,
               mh = 1.67e-24;

  if (ComovingCoordinates != 1) {  
    InitialRedshift = 0; 
    //FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  } 

  /* Quick check to see if sphere overlaps this grid. */

  for (dim = 0; dim < GridRank; dim++)
    if (SphereCenter[dim] - SphereRadius > GridRightEdge[dim] ||
	SphereCenter[dim] + SphereRadius < GridLeftEdge[dim]   )
      return SUCCESS;

  /* Find the units if we are using comoving coordinates. */

  FLOAT DensityConversion = 1, VelocityConversion = 1;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
        TemperatureUnits = 1; 
  FLOAT a = 1, dadt;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in GetUnits.\n");
      return FAIL;
    }

  if (ComovingCoordinates) {  

    CosmologyComputeExpansionFactor(Time, &a, &dadt);
 
    /* Convert cgs units to more reasonable values.
       density:   M(solar)/Mpc^3 
       velocity:  km/s */

    DensityConversion = FLOAT(double(DensityUnits) / SolarMass * pow(Mpc, 3));
    VelocityConversion = FLOAT(double(VelocityUnits) / 1.0e5);

  }
  else {   

    DensityConversion = FLOAT(double(DensityUnits) / SolarMass * pow(Mpc, 3));
    VelocityConversion = FLOAT(double(VelocityUnits) / 1.0e5);

    //DensityConversion = double(DensityUnits);
    //VelocityConversion = double(VelocityUnits);
  }

  /* Compute cell volume and total grid size. */

  FLOAT BoxSize = 1, CellVolume = 1;
  FLOAT DomainWidth[MAX_DIMENSION];

  if (ComovingCoordinates) 
    BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
  else
    BoxSize = LengthUnits/Mpc;  
      
  int size = 1;

  for (dim = 0; dim < GridRank; dim++) {

    CellVolume *= CellWidth[dim][0]*BoxSize;  
    size *= GridDimension[dim];
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }

  FLOAT DConv = CellVolume*DensityConversion;

  /* Set free-free constant, 0.88 is n(e)*n(i) for ionized H/He gas,
     1.15 is a rough value for the mean Gaunt factor, in 10^44 erg/s. */

  FLOAT xray_const1 = 0.88 * 1.15 * 1.43e-27 * pow(DensityUnits / mh, 2) *
                     CellVolume * pow(Mpc, 3) * 1e-44;

  /* This is just CellVolume in cm^3 * 1e-44. The 10^-23 accounts for the
     units used in ComputeXrayEmissivity. */

  FLOAT xray_const2 = CellVolume * pow(Mpc, 3) * 1e-44 * 1e-23;

  /* Find fields: density, total energy, velocity1-3. */

  if (NumberOfBaryonFields > 0) {
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, 
      H2IINum, DINum, DIINum, HDINum;
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
		     HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }

  /* Calculate cooling time */

  float *cooling_time = NULL;
  if (RadiativeCooling) {
    cooling_time = new float[size];
    this->ComputeCoolingTime(cooling_time);
  }

  /* Compute the temperature. */

  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);

  /* If requested, compute emissivity field. */

  float *xray_emissivity = NULL;
  if (parameters->XrayTableFileName != NULL) {
    xray_emissivity = new float[size];
    if (this->ComputeXrayEmissivity(temperature, xray_emissivity,
				    parameters->XrayLowerCutoffkeV,
				    parameters->XrayUpperCutoffkeV,
				    parameters->XrayTableFileName) == FAIL) {
      fprintf(stderr, "error in grid->ComputeXrayEmissivity.\n");
      return FAIL;
    }
  }

  /* Check for metallicty. */

  int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);

  /* Loop over grid. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    if (GridRank > 1) {
      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - SphereCenter[2];
      delz = min(delz, DomainWidth[2]-delz);
    }
    else
      delz = 0;

    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      if (GridRank > 0) {
	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - SphereCenter[1];
	dely = min(dely, DomainWidth[1]-dely);
      }
      else
	dely = 0;
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];

      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
	delx = min(delx, DomainWidth[0]-delx);

	radius2 = delx*delx + dely*dely + delz*delz;

	if (radius2 <= SphereRadius*SphereRadius)   
	  for (n = 0; n < NumberOfBins; n++)
	    if (radius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {  

	      double gas_mass =  
		BaryonField[DensNum][index]*CellVolume*DensityConversion;

	      if (gas_mass == 0)
		break;

	      /* 0) gas density (in Msolar/Mpc^3). */

	      ProfileValue[n][0] += gas_mass; 
	      ProfileWeight[n][0] += CellVolume; 
	      if (ProfileName[0] == NULL) ProfileName[0] = "d_gas (Ms/Mpc^3)";

	      /* 87) clumping factor (gas_density^2). */

	      ProfileValue[n][87] += gas_mass*(gas_mass/CellVolume);
	      ProfileWeight[n][87] += CellVolume;
	      if (ProfileName[87] == NULL) 
		ProfileName[87] = "clumping factor";

	      /* 1) gas rms velocity (in km/s, mass weighted).    
                    (first set vel[] to velocity, correctly adjusted for face-
	             centering if applicable). */

	      for (dim = 0; dim < GridRank; dim++)
		vel[dim] = BaryonField[Vel1Num+dim][index];
	      if (HydroMethod == Zeus_Hydro) {
		vel[0] = 0.5*(vel[0] + BaryonField[Vel1Num][index+1]);
		vel[1] = 0.5*(vel[1] + 
			      BaryonField[Vel2Num][index+GridDimension[0]]);
		vel[2] = 0.5*(vel[2] + 
		BaryonField[Vel3Num][index+GridDimension[0]*GridDimension[1]]);
	      }
	      for (dim = 0; dim < GridRank; dim++)
		vel[dim] = VelocityConversion*vel[dim] - MeanVelocity[dim][0];   

	      for (dim = 0; dim < GridRank; dim++)
		ProfileValue[n][1] += gas_mass * vel[dim] * vel[dim];   
	      ProfileWeight[n][1] += gas_mass;
	      if (ProfileName[1] == NULL) 
		ProfileName[1] = "v_rms_gas_3d (km/s)";

	      /* 3) x-ray luminosity [old:((M(solar)/Mpc^3)^2*sqrt(K))]. 
	            free-free only   [new:  10^44 erg/s] */

#ifdef OLD_WAY
	      float xray = sqrt(temperature[index]) * gas_mass * 
				BaryonField[DensNum][index]*DensityConversion;
	      ProfileWeight[n][3] += CellVolume;
	      if (ProfileName[3] == NULL) 
		ProfileName[3] = "x-ray_lum ((Ms/Mpc^3)^2*K^0.5)";
#else /* OLD_WAY */
	      float xray;
	      if (xray_emissivity == NULL) {
		xray = xray_const1 * BaryonField[DensNum][index] *
		  BaryonField[DensNum][index] * sqrt(temperature[index]);
		ProfileWeight[n][3] = 0;
		if (ProfileName[3] == NULL) 
		  ProfileName[3] = "x-ray_lum (10^44 erg/s)";
	      } else {
		xray = xray_emissivity[index] * xray_const2;
		ProfileWeight[n][3] = 0;
		if (ProfileName[3] == NULL) 
		  ProfileName[3] = "x-ray_luminosity in shell (10^-44 erg/s)";
	      }
#endif /* OLD_WAY */
	      ProfileValue[n][3] += xray;

	      /* 2) gas temperature (in K).  Mass-weighted. */
	      
	      ProfileValue[n][2] += temperature[index] * gas_mass;
	      ProfileWeight[n][2] += gas_mass;
	      if (ProfileName[2] == NULL) 
		ProfileName[2] = "temp_gas_mass_weighted (K)";

	      /* 84) gas temperature (in K).  xray-weighted. */

	      ProfileValue[n][84] += temperature[index] * xray;
	      ProfileWeight[n][84] += xray;
	      if (ProfileName[84] == NULL) 
		ProfileName[84] = "temp_gas_xray_weighted (K)";

	      /* 4) number of samples. */

	      ProfileValue[n][4] += 1.0;
	      ProfileWeight[n][4] = 0;
	      if (ProfileName[4] == NULL) ProfileName[4] = "N_gas";

	      /* 5) mass-weighted entropy. */

	      ProfileValue[n][5] += temperature[index]/
		pow(BaryonField[DensNum][index], float(Gamma-1.0)) * gas_mass;
	      ProfileWeight[n][5] += gas_mass;
	      if (ProfileName[5] == NULL) ProfileName[5] = "S_gas";

	      /* 6) gas radial velocity. */  

	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);
	      ProfileValue[n][6] += radial_vel*gas_mass;
	      ProfileWeight[n][6] += gas_mass;
	      if (ProfileName[6] == NULL) ProfileName[6] = "vr_gas (km/s)";

	      /* 7) gas radial velocity dispersion. */

	      ProfileValue[n][7] += radial_vel*radial_vel*gas_mass;
	      ProfileWeight[n][7] += gas_mass;
	      if (ProfileName[7] == NULL) 
		ProfileName[7] = "vr_rms_gas (km/s)";

	      if (MultiSpecies > 0) {

		/* 8-9) H species. */

		ProfileValue[n][8] += BaryonField[HINum][index]*DConv;
		ProfileWeight[n][8] += gas_mass;
		if (!ProfileName[8]) ProfileName[8] = "HI fraction";
		ProfileValue[n][9] += BaryonField[HIINum][index]*DConv;
		ProfileWeight[n][9] += gas_mass;
		if (!ProfileName[9]) ProfileName[9] = "HII fraction";

		/* 10-12) He species. */

		ProfileValue[n][10] += BaryonField[HeINum][index]*DConv;
		ProfileWeight[n][10] += gas_mass;
		if (!ProfileName[10]) ProfileName[10] = "HeI fraction";
		ProfileValue[n][11] += BaryonField[HeIINum][index]*DConv;
		ProfileWeight[n][11] += gas_mass;
		if (!ProfileName[11]) ProfileName[11] = "HeII fraction";
		ProfileValue[n][12] += BaryonField[HeIIINum][index]*DConv;
		ProfileWeight[n][12] += gas_mass;
		if (!ProfileName[12]) ProfileName[12] = "HeIII fraction";

		/* 13) electron density. */

		ProfileValue[n][13] += BaryonField[HIINum][index]*DConv;
		ProfileWeight[n][13] += gas_mass;
		if (!ProfileName[13]) ProfileName[13] = "e- fraction";

		/* 14-16) H2-related species. */

		if (MultiSpecies > 1) {
		  ProfileValue[n][14] += BaryonField[HMNum][index]*DConv;
		  ProfileWeight[n][14] += gas_mass;
		  if (!ProfileName[14]) ProfileName[14] = "H- fraction";
		  ProfileValue[n][15] += BaryonField[H2INum][index]*DConv;
		  ProfileWeight[n][15] += gas_mass;
		  if (!ProfileName[15]) ProfileName[15] = "H2I fraction";
		  ProfileValue[n][16] += BaryonField[H2IINum][index]*DConv;
		  ProfileWeight[n][16] += gas_mass;
		  if (!ProfileName[16]) ProfileName[16] = "H2II fraction";
		}
		
	      }

	      /* 17-19) angular momentum vector. */  

	      AngWeight = gas_mass * BoxSize;
	      ProfileValue[n][17] += AngWeight * ( vel[1]*delz - vel[2]*dely);  
	      ProfileValue[n][18] += AngWeight * (-vel[0]*delz + vel[2]*delx);
	      ProfileValue[n][19] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	      for (m = 17; m < 20; m++)
		ProfileWeight[n][m] += gas_mass;
	      ProfileName[17] = "Lx_gas (km/s*Mpc)";
	      ProfileName[18] = "Ly_gas (km/s*Mpc)";
	      ProfileName[19] = "Lz_gas (km/s*Mpc)";

	      /* 191-193) Lz_gas/radius -> should be rotational velocity */  

	      AngWeight = gas_mass * BoxSize;
	      ProfileValue[n][191] += AngWeight * ( vel[1]*delz - vel[2]*dely);  
	      ProfileValue[n][192] += AngWeight * (-vel[0]*delz + vel[2]*delx);
	      ProfileValue[n][193] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	      for (m = 191; m < 194; m++)
		ProfileWeight[n][m] += gas_mass * sqrt(radius2);
	      ProfileName[191] = "Lx_gas/radius (km/s)";
	      ProfileName[192] = "Ly_gas/radius (km/s)";
	      ProfileName[193] = "Lz_gas/radius (km/s) = rotational velocity";

	      /* 83) metallicity in gas */

	      if (MetalNum != -1) {
		ProfileValue[n][83] += BaryonField[MetalNum][index]*DConv;
		ProfileWeight[n][83] += gas_mass;
		if (!ProfileName[83]) ProfileName[83] = "Metal Fraction"; 
	      }

	      /* 80-82) angular momentum vector of dense, cold gas. 
	         Only include gas within 0.2 rvir */

	      if (gas_mass/CellVolume > parameters->LowerDensityCutoff &&
		  temperature[index] < parameters->ColdTemperatureCutoff &&
		  radius2 < 0.04*SphereRadius*SphereRadius) {
		if (gas_mass/CellVolume < parameters->UpperDensityCutoff) {
		ProfileValue[n][80] += AngWeight*( vel[1]*delz - vel[2]*dely);
		ProfileValue[n][81] += AngWeight*(-vel[0]*delz + vel[2]*delx);
		ProfileValue[n][82] += AngWeight*( vel[0]*dely - vel[1]*delx);
		}
		for (m = 80; m < 83; m++)
		  ProfileWeight[n][m] += gas_mass;
		ProfileName[80] = "Lx_dense gas (km/s*Mpc)";
		ProfileName[81] = "Ly_dense gas (km/s*Mpc)";
		ProfileName[82] = "Lz_dense gas (km/s*Mpc)";
	      }

	      /* 20) Cooling time */

	      if (RadiativeCooling) {
		ProfileValue[n][20] += cooling_time[index]*TimeUnits * 
		                       gas_mass;
		ProfileWeight[n][20] += gas_mass;
		ProfileName[20] = "T_cool (s)";
	      }

	      /* 120) Compute mu, either from species information or take a
		 guess from temperature (should put in saha...). */

	      if (MultiSpecies > 0) {
		number_density = 
		  0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] + 
			BaryonField[HeIIINum][i]                        ) +
		        BaryonField[HINum][i]   + BaryonField[HIINum][i]  + 
		        BaryonField[DeNum][i];
		if (MultiSpecies > 1)
		  number_density += BaryonField[HMNum][i]   + 
		    0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
		mu = BaryonField[DensNum][i]/number_density;
	      } else {
		if (temperature[index] > 20000)
		  mu = 0.59;
		else
		  mu = 1.22;
	      }
	      ProfileValue[n][120] += mu * gas_mass;
	      ProfileWeight[n][120] += gas_mass;
	      ProfileName[120] = "mu (mh) mean mass per particle";

	      /* 21) Crossing time. */

	      double sound_speed2 = Gamma*kboltz*temperature[index]/
		(mh*mu);
	      
	      //printf ("%d %d %d: T = %g, mh*mu = %g, cs2 = %g, r2 = %g\n", i, j, k,
	      //      temperature[index], mh*mu, sound_speed2, radius2);
	      ProfileValue[n][21] += LengthUnits * sqrt(radius2/sound_speed2) * gas_mass;
	      ProfileWeight[n][21] += gas_mass;
	      ProfileName[21] = "T_cross (s)";

	      /* 23) Cold fraction. */
 
	      if (temperature[index] < parameters->ColdTemperatureCutoff &&
		  gas_mass/CellVolume > parameters->LowerDensityCutoff)
		ProfileValue[n][23] += gas_mass;
	      ProfileWeight[n][23] += gas_mass;
	      ProfileName[23] = "f_cold_dens";
 
	      /* 24) Dense fraction (>DensityCutoff Msolar/Mpc^3). */  
 
	      if (gas_mass/CellVolume > parameters->LowerDensityCutoff)
		ProfileValue[n][24] += gas_mass;
	      ProfileWeight[n][24] += gas_mass;
	      ProfileName[24] = "f_dens";
	      
	      /* 25-27) deuterium-related species. */

	      if (MultiSpecies > 2) {
		ProfileValue[n][25] += BaryonField[DINum][index]*DConv;
		ProfileWeight[n][25] += gas_mass;
		if (!ProfileName[25]) ProfileName[25] = "DI fraction";
		ProfileValue[n][26] += BaryonField[DIINum][index]*DConv;
		ProfileWeight[n][26] += gas_mass;
		if (!ProfileName[26]) ProfileName[26] = "DII fraction";
		ProfileValue[n][27] += BaryonField[HDINum][index]*DConv;
		ProfileWeight[n][27] += gas_mass;
		if (!ProfileName[27]) ProfileName[27] = "HDI fraction";
	      }

	      /* 28) Mass fraction of gas with tcool < thubble */

	      if (RadiativeCooling) {
		if (cooling_time[index] < Time)
		  ProfileValue[n][28] += gas_mass;
		ProfileWeight[n][28] += gas_mass;
		if (!ProfileName[28]) 
		  ProfileName[28] = "f_tcool_lessthan_time";
	      }

	      /* 89) Mass fraction of gas with tcool < thubble and dens>30 
	             Hardcoded to baryon fraction of 10% */

	      if (RadiativeCooling) {
		if (cooling_time[index] < Time && 
		    BaryonField[DensNum][index] > 0.1*30.0)
		  ProfileValue[n][89] += gas_mass;
		ProfileWeight[n][89] += gas_mass;
		if (!ProfileName[89]) 
		  ProfileName[89] = "f_tcool_dens";
		}

	      /* 40-45) inertia tensor. */

	      InertiaWeight = gas_mass * BoxSize * BoxSize;
	      ProfileValue[n][40] += delx * delx * InertiaWeight;
	      ProfileValue[n][41] += delx * dely * InertiaWeight;
	      ProfileValue[n][42] += delx * delz * InertiaWeight;
	      ProfileValue[n][43] += dely * dely * InertiaWeight;
	      ProfileValue[n][44] += dely * delz * InertiaWeight;
	      ProfileValue[n][45] += delz * delz * InertiaWeight;
	      for (m = 40; m < 46; m++)
		ProfileWeight[n][m] += gas_mass;
	      if (ProfileName[40] == NULL) {
		ProfileName[40] = "I_xx_gas";
		ProfileName[41] = "I_xy_gas";
		ProfileName[42] = "I_xz_gas";
		ProfileName[43] = "I_yy_gas";
		ProfileName[44] = "I_yz_gas";
		ProfileName[45] = "I_zz_gas";
	      }

              /* 46) magnitude of vorticity */  

              ip = GRIDINDEX_NOGHOST(i+1,j,k);
              im = GRIDINDEX_NOGHOST(i-1,j,k);
              jp = GRIDINDEX_NOGHOST(i,j+1,k); //#####
              jm = GRIDINDEX_NOGHOST(i,j-1,k);
              kp = GRIDINDEX_NOGHOST(i,j,k+1);
              km = GRIDINDEX_NOGHOST(i,j,k-1);

              ProfileWeight[n][46] += gas_mass;
              if (ProfileName[46] == NULL) ProfileName[46] = "Vorticity (1/s)";

              double vort_x = 1.0, vort_y = 1.0, vort_z = 1.0;

              vort_x = BaryonField[Vel3Num][jp] - BaryonField[Vel3Num][jm] -
                BaryonField[Vel2Num][kp] + BaryonField[Vel2Num][km];
              vort_y = BaryonField[Vel1Num][kp] - BaryonField[Vel1Num][km] -
                BaryonField[Vel3Num][ip] + BaryonField[Vel3Num][im];
              vort_z = BaryonField[Vel2Num][ip] - BaryonField[Vel2Num][im] -
                BaryonField[Vel1Num][jp] + BaryonField[Vel1Num][jm];

              ProfileValue[n][46] += gas_mass * (1.0 / (2*CellWidth[0][i])) *
                sqrt(vort_x*vort_x + vort_y*vort_y + vort_z*vort_z) * TimeUnits;

	      /* Done. */

	      break;  
	    }

      } // end i loop
    } // end j loop
  } // end k loop

  /* Clean up. */

  delete [] temperature;
  if (RadiativeCooling)
    delete [] cooling_time;
  delete [] xray_emissivity;

  } // end: if (NumberOfBaryonFields > 0)

  /* Loop over particles. */

  dely = delz = 0;
  int StarParticle;
 
  for (i = 0; i < NumberOfParticles; i++) {

                      delx = ParticlePosition[0][i] - SphereCenter[0];
    if (GridRank > 1) dely = ParticlePosition[1][i] - SphereCenter[1];
    if (GridRank > 2) delz = ParticlePosition[2][i] - SphereCenter[2];

    radius2 = delx*delx + dely*dely + delz*delz;

    if (radius2 <= SphereRadius*SphereRadius)
      for (n = 0; n < NumberOfBins; n++)
	if (radius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	  float part_mass = ParticleMass[i]*DensityConversion*CellVolume;

	  if (part_mass == 0)
	    break;

	  /* Determine if particle is a star particle. */

	  StarParticle = FALSE;
	  if (NumberOfParticleAttributes > 0 && ParticleType[i] == 2) 
	      StarParticle = TRUE;

	  /* a) Dark matter particles */

	  if (!StarParticle) {

	    /* 30) dm density (in Msolar/Mpc^3). */

	    ProfileValue[n][30] += part_mass;
	    ProfileWeight[n][30] -= CellVolume;
	    if (ProfileName[30] == NULL) ProfileName[30] = "d_dm (Ms/Mpc^3)";

	    /* 31) dm rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][31] += 
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		part_mass;
	    ProfileWeight[n][31] += part_mass;
	    if (ProfileName[31] == NULL) 
	      ProfileName[31] = "v_rms_dm_3d (km/s)";

	    /* 32) number of samples. */

	    ProfileValue[n][32] += 1.0;
	    ProfileWeight[n][32] = 0;
	    if (ProfileName[32] == NULL) ProfileName[32] = "N_dm";

	    /* 33) dm radial velocity. */

	   velx = ParticleVelocity[0][i]*VelocityConversion-MeanVelocity[0][1];
	   vely = ParticleVelocity[1][i]*VelocityConversion-MeanVelocity[1][1];
	   velz = ParticleVelocity[2][i]*VelocityConversion-MeanVelocity[2][1];
	    radial_vel = (delx*velx + dely*vely + delz*velz) / 
	               sqrt(radius2 + 1.0e-20);
	    ProfileValue[n][33] += radial_vel*part_mass;
	    ProfileWeight[n][33] += part_mass;
	    if (ProfileName[33] == NULL) ProfileName[33] = "vr_dm (km/s)";

	    /* 34) dm radial velocity dispersion. */

	    ProfileValue[n][34] += radial_vel*radial_vel*part_mass;
	    ProfileWeight[n][34] += part_mass;
	    if (ProfileName[34] == NULL) ProfileName[34] = "vr_rms_dm (km/s)";

	    /* 35-37) angular momentum vector. */

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][35] += AngWeight * ( vely*delz - velz*dely);
	    ProfileValue[n][36] += AngWeight * (-velx*delz + velz*delx);
	    ProfileValue[n][37] += AngWeight * ( velx*dely - vely*delx);
	    for (m = 35; m < 38; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[35] = "Lx_dm (km/s*Mpc)";
	    ProfileName[36] = "Ly_dm (km/s*Mpc)";
	    ProfileName[37] = "Lz_dm (km/s*Mpc)";

	    /* 50-55) inertia tensor. */

	    InertiaWeight = part_mass * BoxSize * BoxSize;
	    ProfileValue[n][50] += delx * delx * InertiaWeight;
	    ProfileValue[n][51] += delx * dely * InertiaWeight;
	    ProfileValue[n][52] += delx * delz * InertiaWeight;
	    ProfileValue[n][53] += dely * dely * InertiaWeight;
	    ProfileValue[n][54] += dely * delz * InertiaWeight;
	    ProfileValue[n][55] += delz * delz * InertiaWeight;
	    for (m = 50; m < 56; m++)
	      ProfileWeight[n][m] += part_mass;
	    if (ProfileName[50] == NULL) {
	      ProfileName[50] = "I_xx_dm";
	      ProfileName[51] = "I_xy_dm";
	      ProfileName[52] = "I_xz_dm";
	      ProfileName[53] = "I_yy_dm";
	      ProfileName[54] = "I_yz_dm";
	      ProfileName[55] = "I_zz_dm";
	    }

	  } else {

	  /* b) star particles */

	    /* 60) sp density (in Msolar/Mpc^3). */

	    ProfileValue[n][60] += part_mass;
	    ProfileWeight[n][60] -= CellVolume; 
	    if (ProfileName[60] == NULL) ProfileName[60] = "d_sp (Ms/Mpc^3)";

	    /* 61) sp rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][61] += 
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		(ParticleVelocity[dim][i]*VelocityConversion - 
		 MeanVelocity[dim][1])*
		part_mass;
	    ProfileWeight[n][61] += part_mass;
	    if (ProfileName[61] == NULL) 
	      ProfileName[61] = "v_rms_sp_3d (km/s)";

	    /* 62) number of samples. */

	    ProfileValue[n][62] += 1.0;
	    ProfileWeight[n][62] = 0;
	    if (ProfileName[62] == NULL) ProfileName[62] = "N_sp";

	    /* 63) sp radial velocity. */

	   velx = ParticleVelocity[0][i]*VelocityConversion-MeanVelocity[0][1];
	   vely = ParticleVelocity[1][i]*VelocityConversion-MeanVelocity[1][1];
	   velz = ParticleVelocity[2][i]*VelocityConversion-MeanVelocity[2][1];
	    radial_vel = (delx*velx + dely*vely + delz*velz) / 
	               sqrt(radius2 + 1.0e-20);
	    ProfileValue[n][63] += radial_vel*part_mass;
	    ProfileWeight[n][63] += part_mass;
	    if (ProfileName[63] == NULL) ProfileName[63] = "vr_sp (km/s)";

	    /* 64) sp radial velocity dispersion. */

	    ProfileValue[n][64] += radial_vel*radial_vel*part_mass;
	    ProfileWeight[n][64] += part_mass;
	    if (ProfileName[64] == NULL) ProfileName[64] = "vr_rms_sp (km/s)";

	    /* 65-67) angular momentum vector. */

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][65] += AngWeight * ( vely*delz - velz*dely);
	    ProfileValue[n][66] += AngWeight * (-velx*delz + velz*delx);
	    ProfileValue[n][67] += AngWeight * ( velx*dely - vely*delx);
	    for (m = 65; m < 68; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[65] = "Lx_sp (km/s*Mpc)";
	    ProfileName[66] = "Ly_sp (km/s*Mpc)";
	    ProfileName[67] = "Lz_sp (km/s*Mpc)";

	    /* 69) star particle metallicity */

	    ProfileValue[n][69] += part_mass * ParticleAttribute[2][i];
	    ProfileWeight[n][69] += part_mass;
	    if (ProfileName[69] == NULL) ProfileName[69] = "sp_metallicity";

	  } // end if (star particle)

	  /* Done. */

	  break;  

	}  // end of radius loop
  } // end: loop over particles

  return SUCCESS;
}
 
