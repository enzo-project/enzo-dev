/***********************************************************************
/
/  GRID CLASS (ADD THE CONTENTS OF THIS GRID TO THE GIVEN VERTICAL PROFILE)
/
/  written by: Ji-hoon Kim, heavily adopted from Grid_AddToDiskProfile.C
/  date:       December, 2007
/
/  modified1: Ji-hoon Kim
/             July, 2009
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

int grid::AddToVerticalProfile(FLOAT SphereCenter[MAX_DIMENSION], 
			   float SphereRadius,
			   FLOAT MeanVelocity[MAX_DIMENSION][3],
			   int NumberOfBins,
			   FLOAT ProfileRadius[], 
			   FLOAT ProfileValue[][MAX_PROFILES],
			   FLOAT ProfileWeight[][MAX_PROFILES],
			   char *ProfileName[MAX_PROFILES],
			   AnalyzeClusterParameters *parameters,
			   float *DiskVector)
{

  /* Return if the data is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int i, j, k, index, m, n, dim;
  FLOAT radius2, radial_vel, vel[MAX_DIMENSION], circ_vel, AngWeight,
        height, dradius2, length, gas_dens,
        delx, dely, delz, xpos, ypos, zpos;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, mh = 1.67e-24;

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

  float DensityConversion = 1, VelocityConversion = 1;

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

    DensityConversion = float(double(DensityUnits) / SolarMass * pow(Mpc, 3));
    VelocityConversion = float(double(VelocityUnits) / 1.0e5);

  }

  else {   

    DensityConversion = float(double(DensityUnits) / SolarMass * pow(Mpc, 3));
    VelocityConversion = float(double(VelocityUnits) / 1.0e5);

    //    DensityConversion = double(DensityUnits);
    //    VelocityConversion = double(VelocityUnits);

  }

  /* Compute cell volume and total grid size. */

  float BoxSize = 1, CellVolume = 1;
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

  FLOAT DConv = CellVolume*DensityConversion; //Ji-hoon Kim

  int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);  //Ji-hoon Kim

  /* Set free-free constant, 0.88 is n(e)*n(i) for ionized H/He gas,
     1.15 is a rough value for the mean Gaunt factor, in 10^44 erg/s. */

  float xray_const = 0.88 * 1.15 * 1.43e-27 * pow(DensityUnits / mh, 2) *
                     CellVolume * pow(Mpc, 3) * 1e-44;

  /* Compute two vectors which are perpendicular to the disk vector. */

  float theta, phi, axis0[3], axis1[3], tang[3], plane[3], pi = 3.14159;
  theta = atan(sqrt(DiskVector[0]*DiskVector[0] + DiskVector[1]*DiskVector[1])
	       /DiskVector[2]);  // from -pi/2 to +pi/2
  phi   = atan2(DiskVector[1],DiskVector[0]);  // from -pi to +pi

  axis0[0] = sin(theta+0.5*pi)*cos(phi);
  axis0[1] = sin(theta+0.5*pi)*sin(phi);
  axis0[2] = cos(theta+0.5*pi);

  axis1[0] =  (DiskVector[2]*axis0[1] - DiskVector[1]*axis0[2]);
  axis1[1] = -(DiskVector[2]*axis0[0] - DiskVector[0]*axis0[2]);
  axis1[2] =  (DiskVector[1]*axis0[0] - DiskVector[0]*axis0[1]);

#if 0
  printf("DV = %g %g %g\n", DiskVector[0], DiskVector[1], DiskVector[2]);
  printf("a0 = %g %g %g\n", axis0[0], axis0[1], axis0[2]);
  printf("a1 = %g %g %g\n", axis1[0], axis1[1], axis1[2]);
  printf("%g %g %g\n", 
	 DiskVector[0]*axis0[0]+DiskVector[1]*axis0[1]+DiskVector[2]*axis0[2],
	 DiskVector[0]*axis1[0]+DiskVector[1]*axis1[1]+DiskVector[2]*axis1[2],
	 axis1[0]*axis0[0]+axis1[1]*axis0[1]+axis1[2]*axis0[2]);
#endif

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

  /* Compute the temperature. */

  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);

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

	/* If the subgrid field is set, then ignore this cell. */

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;

	float gas_mass =   
	  BaryonField[DensNum][index]*CellVolume*DensityConversion;  

	/* if the cell is not dense and cold, then ignore. */

	if (gas_mass/CellVolume < parameters->LowerDensityCutoff ||       
	    temperature[index] > parameters->ColdTemperatureCutoff)
	  continue;
	
	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - SphereCenter[0];
	delx = min(delx, DomainWidth[0]-delx);

	radius2 = delx*delx + dely*dely + delz*delz;
	if (radius2 > SphereRadius*SphereRadius)
	  continue;

	/* Compute radius in plane. */  

	height = DiskVector[0]*delx + DiskVector[1]*dely + DiskVector[2]*delz;  
	dradius2 = radius2 - height*height;

	if (height*height <= SphereRadius*SphereRadius) 
	  for (n = 0; n < NumberOfBins; n++)
	    if (height*height <= ProfileRadius[n+1]*ProfileRadius[n+1]) { 

	      if (gas_mass == 0)
		break;

	      /* 150) gas density (in Msolar/Mpc^3). */  

	      ProfileValue[n][150] += gas_mass;        
	      ProfileWeight[n][150] += CellVolume;     
	      if (ProfileName[150] == NULL) 
		ProfileName[150] = "d_gas (Ms/Mpc^3)"; 

	      /* 151) gas rms velocity (in km/s, mass weighted). 
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
		ProfileValue[n][151] += gas_mass * vel[dim] * vel[dim];
	      ProfileWeight[n][151] += gas_mass;
	      if (ProfileName[151] == NULL) 
		ProfileName[151] = "v_rms_gas_3d (km/s)";

	      /* 152) gas temperature (in K).  Mass-weighted in cube is
		    less than 40 Mpc (i.e. not a cluster sim), otherwise
	            weight by x-ray luminosity. */

	      ProfileValue[n][152] += temperature[index] * gas_mass;
	      ProfileWeight[n][152] += gas_mass;
	      if (ProfileName[152] == NULL) 
		ProfileName[152] = "temp_gas_mass (K)";
	      
	      /* 153) number of samples. */

	      ProfileValue[n][153] += 1.0;
	      ProfileWeight[n][153] = 0;
	      if (ProfileName[153] == NULL) ProfileName[153] = "N_gas";

	      /* 154) gas radial velocity. */

	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);
	      ProfileValue[n][154] += radial_vel*gas_mass;
	      ProfileWeight[n][154] += gas_mass;
	      if (ProfileName[154] == NULL) ProfileName[154] = "vr_gas (km/s)";

	      /* 155) gas circular velocity. First, compute radial vector 
		 in disk, then compute tangential vector in disk, then
		 normalize it and dot with velocity vector. */

              plane[0] = delx - height*DiskVector[0];  // radial vector in disk
              plane[1] = dely - height*DiskVector[1];
              plane[2] = delz - height*DiskVector[2];
	      tang[0] =  (DiskVector[2]*plane[1] - DiskVector[1]*plane[2]);
	      tang[1] = -(DiskVector[2]*plane[0] - DiskVector[0]*plane[2]);
	      tang[2] =  (DiskVector[1]*plane[0] - DiskVector[0]*plane[1]);
	      length = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]
			    +tiny_number);
	      circ_vel = (tang[0]*vel[0] + tang[1]*vel[1] + tang[2]*vel[2])/
		length;
	      ProfileValue[n][155] += circ_vel*gas_mass;
	      ProfileWeight[n][155] += gas_mass;
	      if (ProfileName[155] == NULL) ProfileName[155] = "vcirc_gas (km/s)";

 	      /* 156) surface density. First just sum the masses, then   
 		 divide by the area at the end */
 
	      ProfileValue[n][156] += gas_mass;
 	      ProfileWeight[n][156] = -1;   
	      if (!ProfileName[156]) ProfileName[156] = "dens_surf (Ms/Mpc^2)";


	      /* ***************************************************************************************
		 160-169 FOR NEWLY ADDED PROFILE FOR VERTICAL PROFILE ONLY (DIFFERENT FROM DISK PROFILE) 
	         *************************************************************************************** */

	      /* 160-162) angular momentum vector. */ 

	      AngWeight = gas_mass * BoxSize; 
	      ProfileValue[n][160] += AngWeight * ( vel[1]*delz - vel[2]*dely);  
	      ProfileValue[n][161] += AngWeight * (-vel[0]*delz + vel[2]*delx);
	      ProfileValue[n][162] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	      for (m = 160; m < 163; m++)
		ProfileWeight[n][m] += gas_mass;
	      ProfileName[160] = "Lx_gas (km/s*Mpc)";
	      ProfileName[161] = "Ly_gas (km/s*Mpc)";
	      ProfileName[162] = "Lz_gas (km/s*Mpc)";

	      /* 163) sum of angular momentum Lz in slab (not mass weighted) */

	      ProfileValue[n][163] += AngWeight * ( vel[0]*dely - vel[1]*delx); //in unit of Ms*km/s*Mpc
 	      ProfileWeight[n][163] = 1;
	      if (!ProfileName[163]) ProfileName[163] = "Total Lz_gas (Ms*km/s*Mpc)";
 
	      /* 164-165) z-directional mass flux 
		          F = integral[ rho(r) * v_z(r) * n_z * dS_xy(r) ] */  

	      /* Below is a quick trick to account for the symmetry up and down the disk.
		 It is very tricky since height vector aligns with DiskVector, which again aligns with L. 
	         If you have a better idea, please change this. */

	      if ( ((vel[2] > 0) && (DiskVector[2] > 0) && (height > 0)) ||  
		   ((vel[2] > 0) && (DiskVector[2] < 0) && (height < 0)) ||
		   ((vel[2] < 0) && (DiskVector[2] > 0) && (height < 0)) ||
		   ((vel[2] < 0) && (DiskVector[2] < 0) && (height > 0)) ) {
		ProfileValue[n][164] += gas_mass/CellVolume 
		                           * abs(vel[2]) * 3.156e7 * 3.241e-20  
		                           * abs(DiskVector[2])  
		                           * CellWidth[0][i] * CellWidth[1][j] * BoxSize * BoxSize;  
		ProfileWeight[n][164] = 1;
		ProfileValue[n][165] += 0;  
		ProfileWeight[n][165] = 1;
	      } else {
		ProfileValue[n][164] += 0;  
		ProfileWeight[n][164] = 1;
		ProfileValue[n][165] += gas_mass/CellVolume 
		                           * abs(vel[2]) * 3.156e7 * 3.241e-20  
		                           * abs(DiskVector[2])  
		                           * CellWidth[0][i] * CellWidth[1][j] * BoxSize * BoxSize;  
		ProfileWeight[n][165] = 1;
	      }	     

	      //velocity in the unit of Mpc/yr, dxdy in (Mpc)^2, Thus: mass flux F in the unit of Ms/yr
 	      if (!ProfileName[164]) ProfileName[164] = "Mass Flux outflow in z direction (Ms/yr)";
	      if (!ProfileName[165]) ProfileName[165] = "Mass Flux  inflow in z direction (Ms/yr)";

	      /* 166) metallicity in gas */

	      if (MetalNum != -1) {
		ProfileValue[n][166] += BaryonField[MetalNum][index]*DConv;
		ProfileWeight[n][166] += gas_mass;
		if (!ProfileName[166]) ProfileName[166] = "Metal Fraction"; 
	      }
 
	      /* Done. */

	      break; 

	    }

      } // end i loop
    } // end j loop
  } // end k loop

  /* Clean up. */

  delete temperature;

  } // end: if (NumberOfBaryonFields > 0)

  /* Loop over particles. */

  dely = delz = 0;
  int StarParticle;
 
  for (i = 0; i < NumberOfParticles; i++) {

                      delx = ParticlePosition[0][i] - SphereCenter[0];
    if (GridRank > 1) dely = ParticlePosition[1][i] - SphereCenter[1];
    if (GridRank > 2) delz = ParticlePosition[2][i] - SphereCenter[2];

    radius2 = delx*delx + dely*dely + delz*delz;
    if (radius2 > SphereRadius*SphereRadius)
      continue;

    /* Compute radius in plane. */

    height = DiskVector[0]*delx + DiskVector[1]*dely + DiskVector[2]*delz;
    dradius2 = radius2 - height*height;


    if (height*height <= SphereRadius*SphereRadius)
      for (n = 0; n < NumberOfBins; n++)
	if (height*height <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	  float part_mass = ParticleMass[i]*DensityConversion*CellVolume;

	  if (part_mass == 0)
	    break;

          /* Determine if particle is a star particle. */

          StarParticle = FALSE;
          if (NumberOfParticleAttributes > 0 && ParticleType[i] == 2) 
              StarParticle = TRUE;

	  /* a) Dark matter particles */

	  if (!StarParticle) {

	    /* 180-182) angular momentum vector. */
	    for (dim = 0; dim < GridRank; dim++)
	      vel[dim] = ParticleVelocity[dim][i]*VelocityConversion - MeanVelocity[dim][1]; 

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][180] += AngWeight * ( vel[1]*delz - vel[2]*dely);
	    ProfileValue[n][181] += AngWeight * (-vel[0]*delz + vel[2]*delx);
	    ProfileValue[n][182] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	    for (m = 180; m < 183; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[180] = "Lx_dm (km/s*Mpc)";
	    ProfileName[181] = "Ly_dm (km/s*Mpc)";
	    ProfileName[182] = "Lz_dm (km/s*Mpc)";

            /* 183) sum of angular momentum Lz in slab (not mass weighted) */
	    ProfileValue[n][183] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	    ProfileWeight[n][183] = 1;
	    if (!ProfileName[183]) ProfileName[183] = "Total Lz_dm (Ms*km/s*Mpc)";

	  } else {

	  /* b) star particles */

	    /* 170) star density (in Msolar/Mpc^3). */

	    ProfileValue[n][170] += part_mass;
	    ProfileWeight[n][170] += CellVolume;  
	    if (ProfileName[170] == NULL) ProfileName[170] = "d_star (Ms/Mpc^3)";

	    /* 171) star rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      vel[dim] = ParticleVelocity[dim][i]*VelocityConversion - MeanVelocity[dim][1]; 

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][171] += vel[dim]*vel[dim]*part_mass;
	    ProfileWeight[n][171] += part_mass;
	    if (ProfileName[171] == NULL) 
	      ProfileName[171] = "v_rms_star_3d (km/s)";

	    /* 172) number of samples. */

	    ProfileValue[n][172] += 1.0;
	    ProfileWeight[n][172] = 0;
	    if (ProfileName[172] == NULL) ProfileName[172] = "N_star";

	    /* 173) stellar circular velocity. First, compute radial vector
                 in disk, then compute tangential vector in disk, then
                 normalize it and dot with velocity vector. */

	    plane[0] = delx - height*DiskVector[0];  // radial vector in disk
	    plane[1] = dely - height*DiskVector[1];
	    plane[2] = delz - height*DiskVector[2];
	    tang[0] =  (DiskVector[2]*plane[1] - DiskVector[1]*plane[2]);
	    tang[1] = -(DiskVector[2]*plane[0] - DiskVector[0]*plane[2]);
	    tang[2] =  (DiskVector[1]*plane[0] - DiskVector[0]*plane[1]);
	    length = sqrt(tang[0]*tang[0]+tang[1]*tang[1]+tang[2]*tang[2]
                            +tiny_number);
	    circ_vel = (tang[0]*vel[0] + tang[1]*vel[1] + tang[2]*vel[2])/
                length;
	    ProfileValue[n][173] += circ_vel*part_mass;
	    ProfileWeight[n][173] += part_mass;
	    if (ProfileName[173] == NULL) ProfileName[173] = "vcirc_star (km/s)";

	    /* 174) stellar surface density (in Msolar/Mpc^2).  */

	    ProfileValue[n][174] += part_mass;
	    ProfileWeight[n][174] = -1;
	    if (ProfileName[174] == NULL) ProfileName[174] = "dens_surf_star (Ms/Mpc^2)";

	    /* 175-177) angular momentum vector. */

	    AngWeight = part_mass * BoxSize;
	    ProfileValue[n][175] += AngWeight * ( vel[1]*delz - vel[2]*dely);
	    ProfileValue[n][176] += AngWeight * (-vel[0]*delz + vel[2]*delx);
	    ProfileValue[n][177] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	    for (m = 175; m < 178; m++)
	      ProfileWeight[n][m] += part_mass;
	    ProfileName[175] = "Lx_sp (km/s*Mpc)";
	    ProfileName[176] = "Ly_sp (km/s*Mpc)";
	    ProfileName[177] = "Lz_sp (km/s*Mpc)";

	    /* 178) sum of angular momentum Lz in slab (not mass weighted) */
	    ProfileValue[n][178] += AngWeight * ( vel[0]*dely - vel[1]*delx);
	    ProfileWeight[n][178] = 1;
	    if (!ProfileName[178]) ProfileName[178] = "Total Lz_sp (Ms*km/s*Mpc)";

	  } // end if (star particle)

	  /* Done. */

	  break; 

	}  // end of radius loop
  } // end: loop over particles

  return SUCCESS;
}
 
