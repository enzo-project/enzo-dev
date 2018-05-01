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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
float gasdev();

static int CollapseTestParticleCount = 0;

int grid::AGNDiskInitializeGrid(float BlackHoleMass,
				int BlackHoleType,
				int DiskType,
				float DiskDensity,
				float DiskTemperature,
				FLOAT DiskRadius,
				FLOAT DiskHeight, 
				int UseGas, int level) 
{
  /* declarations */

  int dim, i, j, k, m, sphere;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* create fields */

  NumberOfBaryonFields = 0;
  int ivel;
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

  if (WritePotential) {
    FieldType[NumberOfBaryonFields++] = GravPotential;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float CosmologySimulationInitialFractionHII   = 1.2e-5;
  float CosmologySimulationInitialFractionHeII  = 1.0e-14;
  float CosmologySimulationInitialFractionHeIII = 1.0e-17;
  float CosmologySimulationInitialFractionHM    = 2.0e-9;
  float CosmologySimulationInitialFractionH2I   = 2.0e-20;
  float CosmologySimulationInitialFractionH2II  = 3.0e-14;

  /* Set various units. */

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.672e-8,
               pi = 3.14159, mh = 1.6726e-24, kboltz = 1.3807e-16;
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  double G = 6.672e-8, Msun = 1.989e33;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits*pow(LengthUnits,3);
  MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;

  /* Set up the baryon field. */
  
  this->AllocateGrids();
   
  /* Loop over the mesh. */

  if (UseGas) {
    
    float density, Velocity[MAX_DIMENSION], temperature, sigma;
    FLOAT x, y = 0, z = 0;
    int n = 0;
    float f_b = 1.0/10.0;
    float RotVelocity[3];
    FLOAT xpos, ypos, cosphi, sinphi, R, Z; 
    double mtot, vrot;
  
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	for (i = 0; i < GridDimension[0]; i++, n++) {
	
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	  density = DiskDensity/1000.0;
	  temperature = DiskTemperature*1000.0;

	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    Velocity[dim] = 0;
	  }
	  
	  /* Find cylindrical coordinates from center. */
	  xpos = x - 0.5;
	  ypos = y - 0.5;
	  
	  R = sqrt(xpos*xpos + ypos*ypos);
	  R = max(R, 0.5*CellWidth[0][0]);
	  Z = z-0.5;
	  cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
	  sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);
	  
	  if (R < DiskRadius && fabs(Z) < 0.5*DiskHeight) { 
	    
	    /* 1) Uniform */
	    
	    if (DiskType == 1) {
	      density = DiskDensity*DiskRadius/(R+DiskHeight);
	      temperature = DiskTemperature;
	      double DiskMass = 2.0*M_PI*DiskDensity*DiskHeight*DiskRadius*(R-DiskHeight*log(R/DiskHeight+1))*MassUnits;
	      double Mass = BlackHoleMass*Msun + DiskMass;
	      printf("BH=%"GSYM", Disk=%"GSYM", r=%"GSYM"\n", BlackHoleMass, DiskMass/Msun, R);
	      vrot = sqrt(G*Mass/(max(R,5*CellWidth[0][0])*LengthUnits))/VelocityUnits;
	      Velocity[0] = -vrot*sinphi;
	      Velocity[1] = vrot*cosphi;
	      
	      // Add thermal noise
	      /*for (dim = 0; dim < GridRank; dim++) {
		Velocity[dim] += gasdev()*sqrt(Gamma*temperature/Mu);
		}*/
	    }

	    /* 2) Circumbinary disk */

	    if (DiskType == 2) {
	      density = DiskDensity;
	      temperature = DiskTemperature;
	      FLOAT a = ExternalGravityRadius/LengthUnits;
	      double Mdisk = DiskDensity*M_PI*R*R*DiskHeight*MassUnits;
	      if (R > 0.5*a) {
		vrot = sqrt(G*(2.0*BlackHoleMass*Msun+Mdisk)/(R*LengthUnits))/VelocityUnits;
	      } else {
		vrot = sqrt(G*(2.0*BlackHoleMass*Msun+Mdisk)/(0.5*a*LengthUnits))*R/(0.5*a)/VelocityUnits;
	      }
	      Velocity[0] = -vrot*sinphi;
	      Velocity[1] = vrot*cosphi;
	    }

	    /* 3) Circumbinary disk with a central gap*/

	    if (DiskType == 3) {
	      FLOAT a = ExternalGravityRadius/LengthUnits;
	      if (R > 1.6*a) {
		density = DiskDensity;
		temperature = DiskTemperature;
	      } else {
		density = DiskDensity/1000.0;
		temperature = DiskTemperature*1000.0;
	      }

	      if (R > 0.5*a) {
		vrot = sqrt(G*2.0*BlackHoleMass*Msun/(R*LengthUnits))/VelocityUnits;
	      } else {
		vrot = sqrt(G*2.0*BlackHoleMass*Msun/(0.5*a*LengthUnits))*R/(0.5*a)/VelocityUnits;
	      }
	      Velocity[0] = -vrot*sinphi;
	      Velocity[1] = vrot*cosphi;
	    }


	    
	  }
	  
	  /* Set density. */
	  
	  BaryonField[iden][n] = density;
	  
	  /* If doing multi-species (HI, etc.), set these. */
	  
	  if (MultiSpecies == 1) {
	    
	    BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
	      CoolData.HydrogenFractionByMass * BaryonField[0][n];
	    BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
	      BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	    BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
	      BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	    BaryonField[HeINum][n] = 
	      (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
	      BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
	  
	    BaryonField[HINum][n] = 
	      CoolData.HydrogenFractionByMass*BaryonField[0][n]
	      - BaryonField[HIINum][n];
	    
	    BaryonField[DeNum][n] = BaryonField[HIINum][n] + 
	    0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	  }
	  
	  /* Set Velocities. */
	  
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[ivel+dim][n] = Velocity[dim];
	  
	  /* Set energy (thermal and then total if necessary). */
	  
	  BaryonField[ietot][n] = temperature/((Gamma-1.0)*Mu);
	  
	  if (DualEnergyFormalism)
	    BaryonField[ieint][n] = BaryonField[ietot][n];
	
	  for (dim = 0; dim < GridRank; dim++) {
	    BaryonField[ietot][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);
	  }
	    
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx][n] = 0.0;
	    BaryonField[iBy][n] = 0.0; 
	    BaryonField[iBz][n] = 556.3/MagneticUnits;
	    BaryonField[iPhi][n] = 0.0;
	    BaryonField[ietot][n] += 0.5*pow(BaryonField[iBz][n],2)/density;
	  }
	  
	} // end for i
      } // end for j
    } // end for k
  } // if (UseGas)

  if (BlackHoleType == 3 && level == MaximumRefinementLevel) {

    double mass_p = BlackHoleMass*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;

    double dxm = dx / pow(2.0, MaximumRefinementLevel);

    NumberOfParticles = 1;
    NumberOfStars = 1;
    //    MaximumParticleNumber = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.5+0.5*dx;
    ParticlePosition[1][0] = 0.5+0.5*dx;
    ParticlePosition[2][0] = 0.5+0.5*dx;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    ParticleAttribute[0][0] = 0.0; // creation time             
    ParticleAttribute[1][0] = 0;
    ParticleAttribute[2][0] = mass_p;
  }

  
  if (BlackHoleType == 1) {

    double mass_p = BlackHoleMass*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;
      

    NumberOfParticles = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_DARK_MATTER;
    ParticlePosition[0][0] = 0.501; // 0.6; // 0.55;                                                    
    ParticlePosition[1][0] = 0.501;
    ParticlePosition[2][0] = 0.501;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    //ParticleAttribute[0][0] = 0.0; // creation time
    //ParticleAttribute[1][0] = t_dyn; // dynamical time
    //ParticleAttribute[2][0] = 0.0; // 

  }

  if (BlackHoleType == 2) {

    double mass_p = 1.0e8*Msun;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);

    NumberOfParticles = 2;
    this->AllocateNewParticles(NumberOfParticles);

    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_DARK_MATTER;
    ParticlePosition[0][0] = 0.55;
    ParticlePosition[1][0] = 0.51;
    ParticlePosition[2][0] = 0.51;

    ParticleMass[1] = den_p;
    ParticleNumber[1] = 1;
    ParticleType[1] = PARTICLE_TYPE_DARK_MATTER;
    ParticlePosition[0][1] = 0.45; 
    ParticlePosition[1][1] = 0.51;
    ParticlePosition[2][1] = 0.51;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = sqrt(0.5*mass_p/0.1);
    ParticleVelocity[2][0] = 0.0;

    ParticleVelocity[0][1] = 0.0;
    ParticleVelocity[1][1] = -sqrt(0.5*mass_p/0.1);
    ParticleVelocity[2][1] = 0.0;

  }


  return SUCCESS;
}
