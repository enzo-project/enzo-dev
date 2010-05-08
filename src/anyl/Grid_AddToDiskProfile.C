/***********************************************************************
/
/  GRID CLASS (ADD THE CONTENTS OF THIS GRID TO THE GIVEN DISK PROFILE)
/
/  written by: Greg Bryan
/  date:       November, 1999
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
void CICDeposit(FLOAT *Image, FLOAT density, float x, float y, int Size);


int grid::AddToDiskProfile(FLOAT SphereCenter[MAX_DIMENSION], 
			   float SphereRadius,
			   FLOAT MeanVelocity[MAX_DIMENSION][3],
			   int NumberOfBins,
			   FLOAT ProfileRadius[], 
			   FLOAT ProfileValue[][MAX_PROFILES],
			   FLOAT ProfileWeight[][MAX_PROFILES],
			   char *ProfileName[MAX_PROFILES],
			   AnalyzeClusterParameters *parameters,
			   float *DiskVector, FLOAT *DiskImage[],
			   int DiskImageSize, float DiskRadius)
{

  /* Return if the data is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int i, j, k, index, n, dim;
  float StarMakerMinimumDynamicalTime, StarMassEjectionFraction;
  double xv1 = 0.0, xv2 = 0.0, minitial = 0.0, mform = 0.0;
  FLOAT dtForSFR;
  FLOAT radius2, radial_vel, vel[MAX_DIMENSION], circ_vel, phi_vel,
        height, dradius2, length, gas_dens,
        delx, dely, delz, xpos, ypos, zpos;
  const double SolarMass = 1.989e33, Mpc = 3.086e24, mh = 1.67e-24, yr = 3.15557e7;

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

  /* Set free-free constant, 0.88 is n(e)*n(i) for ionized H/He gas,
     1.15 is a rough value for the mean Gaunt factor, in 10^44 erg/s. */

  float xray_const = 0.88 * 1.15 * 1.43e-27 * pow(DensityUnits / mh, 2) *
                     CellVolume * pow(Mpc, 3) * 1e-44;

  /* Compute two vectors which are perpindicular to the disk vector. */

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

	if (dradius2 <= SphereRadius*SphereRadius)
	  for (n = 0; n < NumberOfBins; n++)
	    if (dradius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	      if (gas_mass == 0)
		break;

	      /* 100) gas density (in Msolar/Mpc^3). */  

	      ProfileValue[n][100] += gas_mass;
	      ProfileWeight[n][100] += CellVolume;
	      if (ProfileName[100] == NULL) 
		ProfileName[100] = "d_gas (Ms/Mpc^3)";

	      /* 101) gas rms velocity (in km/s, mass weighted). 
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
		ProfileValue[n][101] += gas_mass * vel[dim] * vel[dim];
	      ProfileWeight[n][101] += gas_mass;
	      if (ProfileName[101] == NULL) 
		ProfileName[101] = "v_rms_gas_3d (km/s)";

	      /* 102) gas temperature (in K).  Mass-weighted in cube is
		    less than 40 Mpc (i.e. not a cluster sim), otherwise
	            weight by x-ray luminosity. */

	      ProfileValue[n][102] += temperature[index] * gas_mass;
	      ProfileWeight[n][102] += gas_mass;
	      if (ProfileName[102] == NULL) 
		ProfileName[102] = "temp_gas_mass (K)";
	      
	      /* 103) number of samples. */

	      ProfileValue[n][103] += 1.0;
	      ProfileWeight[n][103] = 0;
	      if (ProfileName[103] == NULL) ProfileName[103] = "N_gas";

	      /* 104) gas radial velocity. */

	      radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);
	      ProfileValue[n][104] += radial_vel*gas_mass;
	      ProfileWeight[n][104] += gas_mass;
	      if (ProfileName[104] == NULL) ProfileName[104] = "vr_gas (km/s)";

	      /* 105) gas circular velocity (v_circ). First, compute radial vector  
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
	      ProfileValue[n][105] += circ_vel*gas_mass;
	      ProfileWeight[n][105] += gas_mass;
	      if (ProfileName[105] == NULL) ProfileName[105] = "vcirc_gas (km/s)";

	      /* 107) gas tangential velocity (v_phi).  JHK in Feb.2008 
		 v_r^2 + v_phi^2 = v^2,  v_phi != v_circ  */
		 
	      phi_vel = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2] - radial_vel*radial_vel);  
	                                                                  
	      ProfileValue[n][107] += phi_vel*gas_mass;              
	      ProfileWeight[n][107] += gas_mass;
	      if (ProfileName[107] == NULL) 
		ProfileName[107] = "vphi_gas (km/s)";

 	      /* 106) surface density. First just sum the masses, then  
 		 divide by the annulus area at the end */
 
	      ProfileValue[n][106] += gas_mass;
 	      ProfileWeight[n][106] = -1;   
	      if (!ProfileName[106]) ProfileName[106] = "dens_surf (Ms/Mpc^2)";
 
	      /* Add to image. */

	      /* project against the two in-disk axis vectors to get
		 index in image. */

	      xpos = delx*axis0[0] + dely*axis0[1] + delz*axis0[2];
	      ypos = delx*axis1[0] + dely*axis1[1] + delz*axis1[2];

	      xpos = (xpos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      ypos = (ypos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      zpos = (height/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
	      
	      gas_dens = gas_mass/
		pow(BoxSize*DiskRadius*SphereRadius/DiskImageSize, 2);

	      if (min(xpos, min(ypos, zpos)) > 0 && 
		  max(xpos, max(ypos, zpos)) < DiskImageSize) {
		CICDeposit(DiskImage[0], gas_dens, xpos, ypos, DiskImageSize);
		CICDeposit(DiskImage[1], gas_dens, xpos, zpos, DiskImageSize);
		CICDeposit(DiskImage[2], gas_dens, ypos, zpos, DiskImageSize);
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

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

 
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



    if (dradius2 <= SphereRadius*SphereRadius)
      for (n = 0; n < NumberOfBins; n++)
	if (dradius2 <= ProfileRadius[n+1]*ProfileRadius[n+1]) {

	  float part_mass = 
	    ParticleMass[i]*DensityConversion*CellVolume;
	  float gas_mass =   
	    BaryonField[DensNum][index]*CellVolume*DensityConversion;  

	  if (part_mass == 0)
	    break;

          /* Determine if particle is a star particle. */

          StarParticle = FALSE;
          if (NumberOfParticleAttributes > 0 && ParticleType[i] == 2) 
              StarParticle = TRUE;

	  /* a) Dark matter particles */

	  if (!StarParticle) {

	  } else {

	  /* b) star particles */

	    /* 110) star density (in Msolar/Mpc^3). */

	    /* The sign of the ProfileWeight below is changed by JHK in Dec.2007
	       For those who wanted to have 'stellar surface density' by this trick, see 114 */

	    ProfileValue[n][110] += part_mass;
	 // ProfileWeight[n][110] -= CellVolume;  
	    ProfileWeight[n][110] += CellVolume;  //JHK

	    if (ProfileName[110] == NULL) ProfileName[110] = "d_star (Ms/Mpc^3)";

	    /* 111) star rms velocity (in km/s, mass weighted). */

	    for (dim = 0; dim < GridRank; dim++)
	      vel[dim] = ParticleVelocity[dim][i]*VelocityConversion - MeanVelocity[dim][1];

	    for (dim = 0; dim < GridRank; dim++)
	      ProfileValue[n][111] += vel[dim]*vel[dim]*part_mass;
	    ProfileWeight[n][111] += part_mass;
	    if (ProfileName[111] == NULL) 
	      ProfileName[111] = "v_rms_star_3d (km/s)";

	    /* 112) number of samples. */

	    ProfileValue[n][112] += 1.0;
	    ProfileWeight[n][112] = 0;
	    if (ProfileName[112] == NULL) ProfileName[112] = "N_star";

	    /* 113) stellar circular velocity. First, compute radial vector
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
	    ProfileValue[n][113] += circ_vel*part_mass;
	    ProfileWeight[n][113] += part_mass;
	    if (ProfileName[113] == NULL) ProfileName[113] = "vcirc_star (km/s)";

	    /* 115) stellar tangential velocity.  JHK in Feb.2008 
	       see 'gas tangential velocity' for the definition */

	    radial_vel = (delx*vel[0] + dely*vel[1] + delz*vel[2]) / 
		           sqrt((radius2 == 0)? 1.0 : radius2);
	    phi_vel = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2] - radial_vel*radial_vel);  
	                                                                
	    ProfileValue[n][115] += phi_vel*part_mass;              
	    ProfileWeight[n][115] += part_mass;
	    if (ProfileName[115] == NULL) ProfileName[115] = "vphi_star (km/s)";

	    /* 114) stellar surface density (in Msolar/Mpc^2) as in 106.  
	       JHK in Dec.2007 */

	    ProfileValue[n][114] += part_mass;
	    ProfileWeight[n][114] = -1;
	    if (ProfileName[114] == NULL) ProfileName[114] = "dens_surf_star (Ms/Mpc^2)";

	    /* 116) SFR surface density. First just sum the forming stellar masses, 
	       then divide by the annulus area at the end.  This assumes stars were 
	       created with Cen & Ostriker algorithm (check star_maker7).  JHK in Jan.2010 */  

	    StarMakerMinimumDynamicalTime = 1e6; //yr
	    StarMassEjectionFraction      = 0.8; //#####

	    dtForSFR = StarMakerMinimumDynamicalTime * yr / TimeUnits;
	    xv1 = (Time            - ParticleAttribute[0][i])/ParticleAttribute[1][i];
	    xv2 = (Time + dtForSFR - ParticleAttribute[0][i])/ParticleAttribute[1][i];

	    /*
	    if (xv1 < 1) {
	      ProfileValue[n][116] += part_mass / ( dtForSFR * TimeUnits / yr );
	      ProfileWeight[n][116] = -1;   
	      if (ProfileName[116] == NULL) ProfileName[116] = "dens_surf_SFR (Ms/yr/Mpc^2)";
	    }	
	    */

	    /* if the star is not active or if the gas is not dense, ignore. */

	    if (xv2 < 12 && gas_mass/CellVolume > parameters->LowerDensityCutoff) {
	      minitial = part_mass / (1.0 - StarMassEjectionFraction*(1.0 - (1.0 + xv1)*exp(-xv1)));
	      mform = minitial * ((1.0 + xv1)*exp(-xv1) - (1.0 + xv2)*exp(-xv2));
	      mform = max(min(mform, part_mass), 0.0);

	      /*
	      if (i < 10000) {
	      printf("ParticleAttribute[0][i] = %g, ParticleAttribute[1][i] = %g, dtForSFR = %g\n", 
		     ParticleAttribute[0][i], ParticleAttribute[1][i], dtForSFR);
	      printf("xv1 = %g, xv2 = %g, part_mass = %g, minitial = %g, mform = %g\n", 
		     xv1, xv2, part_mass, minitial, mform);
	      } 
	      */

	      ProfileValue[n][116] += (1.0 - StarMassEjectionFraction) * mform / ( dtForSFR * TimeUnits / yr );
	      ProfileWeight[n][116] = -1; 
	      if (ProfileName[116] == NULL) ProfileName[116] = "dens_surf_SFR (Ms/yr/Mpc^2)";
	    }	      

	    /* Add to image. */

	    /* project against the two in-disk axis vectors to get
	       index in image. */

              xpos = delx*axis0[0] + dely*axis0[1] + delz*axis0[2];
              ypos = delx*axis1[0] + dely*axis1[1] + delz*axis1[2];

              xpos = (xpos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
              ypos = (ypos/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;
              zpos = (height/(DiskRadius*SphereRadius)+0.5)*DiskImageSize;

              gas_dens = part_mass/
                pow(BoxSize*DiskRadius*SphereRadius/DiskImageSize, 2);

              if (min(xpos, min(ypos, zpos)) > 0 &&
                  max(xpos, max(ypos, zpos)) < DiskImageSize) {
                CICDeposit(DiskImage[3], gas_dens, xpos, ypos, DiskImageSize);
                CICDeposit(DiskImage[4], gas_dens, xpos, zpos, DiskImageSize);
                CICDeposit(DiskImage[5], gas_dens, ypos, zpos, DiskImageSize);
              }

	  } // end if (star particle)

	  /* Done. */

	  break; 
	}  // end of radius loop
  } // end: loop over particles

  return SUCCESS;
}
 
void CICDeposit(FLOAT *Image, FLOAT density, float x, float y, int Size)
{
  int i, j, ip1, jp1;
  float dx, dy;

  i = int(x + 0.5) - 1;
  j = int(y + 0.5) - 1;

  if (i < 0 || j < 0 || i >= Size || j >= Size)
    return;

  dx = i + 1.5 - x;
  dy = j + 1.5 - y;

  ip1 = min(i+1, Size-1);
  jp1 = min(j+1, Size-1);

  Image[i  +j  *Size] += density*     dx *     dy;
  Image[ip1+j  *Size] += density*(1.0-dx)*     dy;
  Image[i  +jp1*Size] += density*     dx *(1.0-dy);
  Image[ip1+jp1*Size] += density*(1.0-dx)*(1.0-dy);

}
