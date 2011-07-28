/***********************************************************************
/
/  CALCULATE MASS ACCRETION ONTO STAR PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: Ji-hoon Kim
/             July, 2009
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::CalculateMassAccretion(float &BondiRadius, float &density)
{

  if ((this->type != BlackHole && this->type != MBH) || 
      (this->CurrentGrid == NULL)) 
    return SUCCESS;

  if (this->type == MBH && MBHAccretion == 0)
    return SUCCESS;

  const double Grav = 6.673e-8, k_b = 1.38e-16, m_h = 1.673e-24, Mpc = 3.086e24;
  const double Msun = 1.989e33, yr = 3.1557e7, sigma_T = 6.65e-25, c = 3.0e10;
  const int AccretionType = LOCAL_ACCRETION;
  FLOAT time = CurrentGrid->OldTime;

  if (time <= 0)
    time = CurrentGrid->Time - CurrentGrid->dtFixed;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					      Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (CurrentGrid->
	IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) 
	== FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }

  int igrid[MAX_DIMENSION], dim, index, size = 1;
  float c_s, mu, number_density, old_mass, delta_mass, mdot;
  float mdot_original = 0.0, mdot_UpperLimit, mdot_Edd;
  float *temperature, v_rel, dvel;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= CurrentGrid->GridDimension[dim];
    igrid[dim] = (int) ((pos[dim] - CurrentGrid->GridLeftEdge[dim]) /
			CurrentGrid->CellWidth[0][0]);
  }

  temperature = new float[size];
  if (CurrentGrid->ComputeTemperatureField(temperature) == FAIL) {
    ENZO_FAIL("Error in ComputeTemperatureField.\n");
  }

  /* Reset the accretion rate (DeltaMass) */

  this->ResetAccretion();

  if (AccretionType == LOCAL_ACCRETION && this->type != MBH) {

    // Calculate gas density inside cell
    index = 
      ((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
       igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
      igrid[0] + CurrentGrid->GridStartIndex[0];
    density = CurrentGrid->BaryonField[DensNum][index];
    if (MultiSpecies == 0) {
      number_density = density * DensityUnits / (Mu * m_h);
      mu = Mu;
    } else {
      number_density = 
	CurrentGrid->BaryonField[HINum][index] + 
	CurrentGrid->BaryonField[HIINum][index] +
	CurrentGrid->BaryonField[DeNum][index] +
	0.25 * (CurrentGrid->BaryonField[HeINum][index] +
		CurrentGrid->BaryonField[HeIINum][index] +
		CurrentGrid->BaryonField[HeIIINum][index]);
      if (MultiSpecies > 1)
	number_density += 
	  CurrentGrid->BaryonField[HMNum][index] +
	  0.5 * (CurrentGrid->BaryonField[H2INum][index] +
		 CurrentGrid->BaryonField[H2IINum][index]);
      mu = density / number_density;
    }

    c_s = sqrt(Gamma * k_b * temperature[index] / (mu * m_h));  
    old_mass = (float)(this->Mass);

    // Calculate gas relative velocity (cm/s)
    v_rel = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = vel[dim] - CurrentGrid->BaryonField[Vel1Num+dim][index];
      v_rel += delta_vel[dim] * delta_vel[dim];
    }
    v_rel = sqrt(v_rel) * VelocityUnits;

    // Calculate accretion rate in Msun/s
    mdot = 4.0 * PI * Grav*Grav * (old_mass * old_mass * Msun) * 
      (density * DensityUnits) / pow(c_s * c_s + v_rel * v_rel, 1.5);

    //this->DeltaMass += mdot * (CurrentGrid->dtFixed * TimeUnits);

    this->naccretions = 1;
    if (this->accretion_rate == NULL) {
      this->accretion_rate = new float[naccretions];
      this->accretion_time = new FLOAT[naccretions];
    }
    this->accretion_rate[0] = mdot;
    this->accretion_time[0] = time;

    if (mdot > 0.0)
      fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
	      "M_BH = %lf Msun, rho = %"GSYM" g/cm3, T = %"GSYM" K, v_rel = %"GSYM" cm/s\n",
	      Identifier, time, mdot*yr, Mass, density*DensityUnits,
	      temperature[index], v_rel);

  } // ENDIF LOCAL_ACCRETION with type != MBH






  /* Let us just make a separate routine for MBH Accretion;
     It became too comlicated...  By Ji-hoon Kim in Dec.2010 */

  if (AccretionType == LOCAL_ACCRETION && this->type == MBH) {

    // Calculate gas density inside cell
    index = 
      ((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
       igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
      igrid[0] + CurrentGrid->GridStartIndex[0];
    density = CurrentGrid->BaryonField[DensNum][index];

//     fprintf(stdout, "index = %d, density = %g\n", index, density);  
//     fprintf(stdout, "igrid[0], igrid[1], igrid[2] = %d, %d, %d\n", igrid[0], igrid[1], igrid[2]); 
//     fprintf(stdout, "pos[0], [1], [2] = %g %g %g\n", pos[0], pos[1], pos[2]); 
//     fprintf(stdout, "GridLeftEdge[0], [1], [2] = %g, %g, %g\n",
// 	    CurrentGrid->GridLeftEdge[0], CurrentGrid->GridLeftEdge[1], CurrentGrid->GridLeftEdge[2]);
//     fprintf(stdout, "GridDimension[0], [1], [2] = %d, %d, %d\n",
// 	    CurrentGrid->GridDimension[0], CurrentGrid->GridDimension[1], CurrentGrid->GridDimension[2]);

    if (MultiSpecies == 0) {
      number_density = density * DensityUnits / (Mu * m_h);
      mu = Mu;
    } else {
      /*
      printf("star::CMA: HI, HII, De, HeI, HeII, HeIII, index = %d %d %d %d %d %d %d\n",
	     HINum, HIINum, DeNum, HeINum, HeIINum, HeIIINum, index);
      */
      number_density = 
	CurrentGrid->BaryonField[HINum][index] + 
	CurrentGrid->BaryonField[HIINum][index] +
	CurrentGrid->BaryonField[DeNum][index] +
	0.25 * (CurrentGrid->BaryonField[HeINum][index] +
		CurrentGrid->BaryonField[HeIINum][index] +
		CurrentGrid->BaryonField[HeIIINum][index]);
      if (MultiSpecies > 1)
	number_density += 
	  CurrentGrid->BaryonField[HMNum][index] +
	  0.5 * (CurrentGrid->BaryonField[H2INum][index] +
		 CurrentGrid->BaryonField[H2IINum][index]);
      mu = density / number_density;
    }

    // Calculate gas relative velocity (cm/s)
    v_rel = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      delta_vel[dim] = vel[dim] - CurrentGrid->BaryonField[Vel1Num+dim][index];
      v_rel += delta_vel[dim] * delta_vel[dim];
    }
    v_rel = sqrt(v_rel) * VelocityUnits;
 
    // For MBH, if requested, fix the temperature and/or zero v_rel 
    // so you don't get overpowered by high T SN bubble

    if (MBHAccretion == 2 || MBHAccretion == 12) {
      v_rel = 0.0;
      temperature[index] = MBHAccretionFixedTemperature;   
    } else if (MBHAccretion == 4 || MBHAccretion == 14 ||
	       MBHAccretion == 5 || MBHAccretion == 15) {
      v_rel = 0.0;
    }

    c_s = sqrt(Gamma * k_b * temperature[index] / (mu * m_h));
    old_mass = (float)(this->Mass);

    // Calculate accretion rate in Msun/s
    mdot = 4.0 * PI * Grav*Grav * (old_mass * old_mass * Msun) * 
      (density * DensityUnits) / pow(c_s * c_s + v_rel * v_rel, 1.5);

    /* For MBH, MBHAccretingMassRatio is implemented; when we resolve Bondi radius, 
       the local density used to calculate mdot can be higher than what was supposed 
       to be used (density at Bondi radius), resulting in the overestimation of mdot. 
       MBHAccretingMassRatio (=< 1) can be used to fix this.  Or one might try using
       the density profile of R^-1.5 to estimate the density at R_Bondi 
       (similar to Wang et al. 2009) -Ji-hoon Kim, Sep.2009 */

    mdot_original = mdot;
    BondiRadius = 2.0 * Grav * old_mass * Msun / (c_s * c_s + v_rel * v_rel) / LengthUnits;

    if (MBHAccretingMassRatio > 0) {
      mdot *= MBHAccretingMassRatio;
    } else if (MBHAccretingMassRatio == BONDI_ACCRETION_CORRECT_ANALYTIC) {
      mdot *= min(pow(CurrentGrid->CellWidth[0][0]/BondiRadius, 1.5), 1.0);
    } 
    
    /* Don't take out too much mass suddenly; this is usually not needed 
       because mdot is almost always very small */

//  mdot_UpperLimit = 0.10 * density * DensityUnits * 
//	pow(CurrentGrid->CellWidth[0][0]*LengthUnits, 3.0) / Msun / 
//	(CurrentGrid->dtFixed) / TimeUnits;
//  mdot = min(mdot, mdot_UpperLimit);

      
    /* If requested, just fix mdot (e.g. to 1e-4 Msun/yr) */

    if (MBHAccretion == 3 || MBHAccretion == 13)
      mdot = MBHAccretionFixedRate / yr; 

    /* If requested, modify (suppress) Bondi-Hoyle formalism 
       by taking vorticity into account using Krumholz et al.(2006); see Eq.(3) */

    if (MBHAccretion == 5 || MBHAccretion == 15) {

      if (CurrentGrid->GridRank != 3) 
	ENZO_FAIL("Devised only for three dimension.");

      double vorticity;
      float curl_x, curl_y, curl_z;
      int index_yp1, index_ym1, index_zp1, index_zm1;
      FLOAT dx = CurrentGrid->CellWidth[0][0], 
	dy = CurrentGrid->CellWidth[1][0], 
	dz = CurrentGrid->CellWidth[2][0];

      index_yp1 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + 1 + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_ym1 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] - 1 + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_zp1 = 
	((igrid[2] + 1 + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_zm1 = 
	((igrid[2] - 1 + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];      
      curl_x = 
	float((0.5*(CurrentGrid->BaryonField[Vel3Num][index_yp1] - CurrentGrid->BaryonField[Vel3Num][index_ym1])/dy -
	       0.5*(CurrentGrid->BaryonField[Vel2Num][index_zp1] - CurrentGrid->BaryonField[Vel2Num][index_zm1])/dz));
      curl_y = 
	float((0.5*(CurrentGrid->BaryonField[Vel1Num][index_zp1] - CurrentGrid->BaryonField[Vel1Num][index_zm1])/dz -
	       0.5*(CurrentGrid->BaryonField[Vel3Num][index+1] - CurrentGrid->BaryonField[Vel3Num][index-1])/dx));
      curl_z = 
	float((0.5*(CurrentGrid->BaryonField[Vel2Num][index+1] - CurrentGrid->BaryonField[Vel2Num][index-1])/dx -
	       0.5*(CurrentGrid->BaryonField[Vel1Num][index_yp1] - CurrentGrid->BaryonField[Vel1Num][index_ym1])/dy));

      // Calculate dimensionless vorticity described in Krumholz et al.
      vorticity = sqrt(curl_x*curl_x + curl_y*curl_y + curl_z*curl_z) / TimeUnits * 
	BondiRadius * LengthUnits / c_s;

      mdot *= pow(1 + pow(0.34 / (1 + pow(vorticity, 0.9)), -2.0), -0.5);

      if (mdot > 0.0) {
	fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
		"M_BH = %lf Msun, rho = %"GSYM" g/cm3, c_s = %"GSYM" cm/s, "
		"vorticity = %"GSYM" /s, suppression factor = %"GSYM"\n",
		Identifier, time, mdot*yr, Mass, density*DensityUnits, c_s, vorticity,
		pow(1 + pow(0.34 / (1 + pow(vorticity, 0.9)), -2.0), -0.5));
//	this->PrintInfo();  
      }

    }

    /* If requested, adandon Bondi-Hoyle formalism and
       use the alpha disk formalism in DeBuhr et al.(2010); see Eq.(6)  */

    if (MBHAccretion == 6 || MBHAccretion == 16) {

      double alpha = 0.1; 

      // Calculate accretion rate in Msun/s
      // mdot = 3pi*alpha * c_s^2 * Sigma / Omega
      //      = 3pi*alpha * c_s^2 * rho*R / sqrt(G*M/(R/2)^3)
      //      = 3pi*alpha * c_s^2 * rho*R / sqrt(G*rho*8 + G*M_BH/(R/2)^3)
      mdot = 3.0 * PI * alpha * (c_s * c_s) / Msun * 
	(density * DensityUnits * CurrentGrid->CellWidth[0][0]*LengthUnits) /
	sqrt(Grav * density * DensityUnits * 8.0 + 
	     Grav * (old_mass * Msun) / pow(CurrentGrid->CellWidth[0][0]*LengthUnits/2.0, 3.0));

      if (mdot > 0.0) {
	fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
		"M_BH = %lf Msun, rho = %"GSYM" g/cm3, c_s = %"GSYM" cm/s, T = %"GSYM" K, "
                "Omega1 = %"GSYM" /s, Omeag2 = %"GSYM" /s\n",
		Identifier, time, mdot*yr, Mass, density*DensityUnits, c_s, temperature[index],
		sqrt(Grav * density * DensityUnits * 8.0), 	     
		sqrt(Grav * old_mass * Msun / pow(CurrentGrid->CellWidth[0][0]*LengthUnits/2.0, 3.0)));
//	this->PrintInfo();  
      }

    }		   

    /* If requested, adandon Bondi-Hoyle formalism and
       use modified version of the alpha disk formalism derived with Cen (2011) */

    if (MBHAccretion == 7 || MBHAccretion == 17) {

      int index_L;
      double alpha  = 0.1, kappa  = 0.4, lambda;
      double gas_angmom[] = {0.0, 0.0, 0.0}, total_gas_mass = 0.0, gas_mass = 0.0;
      FLOAT CellVolume = 1, BoxSize = 1, DensityConversion = 1, VelocityConversion = 1;
      FLOAT a = 1, dadt;
      FLOAT delx, dely, delz, velx, vely, velz;
      const double sigma_SB = 5.67e-5;

      if (ComovingCoordinates) {
	CosmologyComputeExpansionFactor(time, &a, &dadt);
	BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
      } else
	BoxSize = LengthUnits/Mpc; // to Mpc

      DensityConversion = FLOAT(double(DensityUnits) / Msun * pow(Mpc, 3)); // to Msun/Mpc^3
      VelocityConversion = FLOAT(double(VelocityUnits) / 1.0e5); // to km/s

      for (dim = 0; dim < MAX_DIMENSION; dim++) 
	CellVolume *= CurrentGrid->CellWidth[0][0]*BoxSize; // in Mpc^3

      /* Find angular momentum in 27 cells */
	
      for (int kk = -1; kk <= 1; kk++) {
	// relative position
	delz = (CurrentGrid->CellLeftEdge[2][igrid[2] + kk + CurrentGrid->GridStartIndex[2]] + 
		0.5*CurrentGrid->CellWidth[2][0] - pos[2]) * BoxSize; // in Mpc
	
	for (int jj = -1; jj <= 1; jj++) {
	  dely = (CurrentGrid->CellLeftEdge[1][igrid[1] + jj + CurrentGrid->GridStartIndex[1]] + 
		  0.5*CurrentGrid->CellWidth[1][0] - pos[1]) * BoxSize;
	  
	  for (int ii = -1; ii <= 1; ii++) {
	    delx = (CurrentGrid->CellLeftEdge[0][igrid[0] + ii + CurrentGrid->GridStartIndex[0]] + 
		    0.5*CurrentGrid->CellWidth[0][0] - pos[0]) * BoxSize;

	    index_L = 
	      ((igrid[2] + kk + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	       igrid[1] + jj + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	      igrid[0] + ii + CurrentGrid->GridStartIndex[0];

	    gas_mass = CurrentGrid->BaryonField[DensNum][index_L] * 
	      DensityConversion * CellVolume; // in Msun
	    
	    // relative velocity
	    velx = (CurrentGrid->BaryonField[Vel1Num][index_L] - vel[0]) * VelocityConversion; // in km/s
	    vely = (CurrentGrid->BaryonField[Vel2Num][index_L] - vel[1]) * VelocityConversion;
	    velz = (CurrentGrid->BaryonField[Vel3Num][index_L] - vel[2]) * VelocityConversion;
	    
	    // store gas angular momentum in: Msun * Mpc * km/s
	    gas_angmom[0] += gas_mass * ( vely*delz - velz*dely); 
	    gas_angmom[1] += gas_mass * (-velx*delz + velz*delx);
	    gas_angmom[2] += gas_mass * ( velx*dely - vely*delx);
	    total_gas_mass += gas_mass;
	  }
	}
      }

      // specific gas angular momentum in: Mpc * km/s
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	gas_angmom[dim] /= total_gas_mass;

      // the specific angular momentum in Keplerian orbit in: Mpc * km/s
      double Keplerian_angmom = sqrt(Grav * old_mass * Msun / (CurrentGrid->CellWidth[0][0]*LengthUnits)) / 1e5 *
	CurrentGrid->CellWidth[0][0] * BoxSize;

      // now lambda = the ratio of specific angular momentum (the real vs. the Keplerian orbit)
      lambda = fabs(gas_angmom[2]) / Keplerian_angmom;

      if (mdot > 0.0) {
 	fprintf(stdout, "Star::CalculateMassAccretion: method=7: lambda = %g, gas_angmom[2] = %g, "
 		"Keplerian_angmom = %g\n", lambda, gas_angmom[2], Keplerian_angmom);
      }



#ifdef DONT_BOTHER_AND_JUST_FIX_LAMBDA
      lambda = 0.1;
#endif

      // Calculate accretion rate in Msun/s
      // mdot = 3pi*alpha * c_s^2 * Sigma / Omega
      //      = 3^(4/3)/4^(1/3) * pi * (3*k_b/m_h)^(4/3) * (3*kappa/(16*sigma_SB*G))^(1/3) * alpha^(4/3) *
      //        M^(-1/3) * Sigma^(5/3) * R * lambda^(-4/3)
      mdot = pow(3.0, 4.0/3.0)/pow(4.0, 1.0/3.0) * PI / Msun * pow(3.0*k_b/m_h, 4.0/3.0) *
	pow(3.0*kappa/16.0/sigma_SB/Grav, 1.0/3.0) * pow(alpha, 4.0/3.0) * pow(old_mass * Msun, -1.0/3.0) *
	pow(density * DensityUnits * CurrentGrid->CellWidth[0][0]*LengthUnits, 5.0/3.0) *
	CurrentGrid->CellWidth[0][0]*LengthUnits * pow(lambda, -7.0/3.0);

//    if (mdot > 0.0) {
// 	fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
// 		"M_BH = %lf Msun, rho = %"GSYM" g/cm3, c_s = %"GSYM" cm/s, T = %"GSYM" K\n",
// 		Identifier, time, mdot*yr, Mass, density*DensityUnits, c_s, temperature[index]);
//	this->PrintInfo();  
//    }

    }		   

    /* If requested, adandon Bondi-Hoyle formalism and
       just use the actual converging mass into the cell to calculate mdot */

    if (MBHAccretion == 8 || MBHAccretion == 18) {

      int index_i1, index_j1, index_k1, index_i2, index_j2, index_k2, divergence;
      index_i1 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + 1 + CurrentGrid->GridStartIndex[0];
      index_j1 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + 1 + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_k1 = 
	((igrid[2] + 1 + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_i2 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] - 1 + CurrentGrid->GridStartIndex[0];
      index_j2 = 
	((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] - 1 + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];
      index_k2 = 
	((igrid[2] - 1 + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
	 igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
	igrid[0] + CurrentGrid->GridStartIndex[0];

      if (HydroMethod == 2) {
	divergence = CurrentGrid->BaryonField[Vel1Num+dim][index_i1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index]
	  + CurrentGrid->BaryonField[Vel1Num+dim][index_j1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index]
	  + CurrentGrid->BaryonField[Vel1Num+dim][index_k1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index];
      } else {
	divergence = CurrentGrid->BaryonField[Vel1Num+dim][index_i1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index_i2]
	  + CurrentGrid->BaryonField[Vel1Num+dim][index_j1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index_j2]
	  + CurrentGrid->BaryonField[Vel1Num+dim][index_k1]
	  - CurrentGrid->BaryonField[Vel1Num+dim][index_k2];
      }

      // Calculate accretion rate in Msun/s
      // mdot = -rho * A * div(v) in Msun/sec (if no converging flow, no accretion)
      mdot = max(0.0, -density * DensityUnits * pow(CurrentGrid->CellWidth[0][0]*LengthUnits, 2.0) *
		 divergence * VelocityUnits / Msun);	
    }		   

    /* Calculate Eddington accretion rate in Msun/s; In general, when calculating 
       the Eddington limit the radiative efficiency shouldn't be smaller than 
       0.1 even when MBHFeedbackRadiativeEfficiency is set to be lower than 0.1 for 
       a logistical purpose (to combine radiative+mechanical feedbacks for example) */

    mdot_Edd = 4.0 * PI * Grav * old_mass * m_h /
      max(MBHFeedbackRadiativeEfficiency, 0.1) / sigma_T / c;     
    if (MBHAccretion < 10 && MBHAccretingMassRatio != BONDI_ACCRETION_CORRECT_NUMERICAL) {
      mdot = min(mdot, mdot_Edd); 
    }

    /* No accretion if the BH is in some low-density and cold cell. */
    if (density < tiny_number || temperature[index] < 10 || isnan(mdot)) 
      mdot = 0.0;

    //this->DeltaMass += mdot * (CurrentGrid->dtFixed * TimeUnits);
    
    this->naccretions = 1;
    if (this->accretion_rate == NULL) {
      this->accretion_rate = new float[naccretions];
      this->accretion_time = new FLOAT[naccretions];
    }
    this->accretion_rate[0] = mdot;
    this->accretion_time[0] = time;
    
    if (mdot > 0.0) {
      fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" (%"GSYM"/%"GSYM") Msun/yr, "
	      "M_BH = %lf Msun, rho = %"GSYM" g/cm3, T = %"GSYM" K, v_rel = %"GSYM" cm/s\n",
	      Identifier, time, mdot*yr, mdot_original*yr, mdot_Edd*yr, Mass, density*DensityUnits,
	      temperature[index], v_rel);
      //    this->PrintInfo();  
    }
    
  } // ENDIF LOCAL_ACCRETION with type == MBH







  if (AccretionType == BONDI_ACCRETION ||
      AccretionType == RADIAL_ACCRETION) {
    ENZO_VFAIL("AccretionType = %"ISYM" not implemented yet.\n", AccretionType);
  }

  delete [] temperature;

  return SUCCESS;
}
