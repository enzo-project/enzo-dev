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
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double Gaussian(double cs);
double BE(double r);

float InitialFractionHII = 0.0;
float InitialFractionHeII = 0.0;
float InitialFractionHeIII = 0.0;

int grid::Collapse3DInitializeGrid(int n_sphere,
				   FLOAT r_sphere[MAX_SPHERES],
				   FLOAT rc_sphere[MAX_SPHERES],
				   float rho_sphere[MAX_SPHERES],
				   float p_sphere[MAX_SPHERES],
				   float cs_sphere[MAX_SPHERES],
				   FLOAT sphere_position[MAX_SPHERES][MAX_DIMENSION],
				   float omega_sphere[MAX_SPHERES],
				   int   sphere_type[MAX_SPHERES],
				   float rho_medium, float p_medium, int level)
{
  /* declarations */

  int dim, i, j, k, m, sphere;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum, 
    kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3;


  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
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


  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }
#ifdef TRANSFER
  if (RadiativeTransfer && (MultiSpecies < 1)) {
    fprintf(stderr, "Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
    return FAIL;
  }
  //   Allocate fields for photo ionization and heating rates                           
  if (RadiativeTransfer)
    if (MultiSpecies) {
      FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
      FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
      FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
      if (MultiSpecies > 1) {
        FieldType[kdissH2INum    = NumberOfBaryonFields++] = kdissH2I;
      }
    }

  if (RadiationPressure && RadiativeTransfer) {
    FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
    FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
    FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
  }
  NumberOfPhotonPackages = 0;
  PhotonPackages-> NextPackage= NULL;
#endif


  float rhou, lenu, tempu, tu, velu;  
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
  double massu = rhou*pow(lenu,3);


  this->AllocateGrids();

  /* Initialize radiation fields */
#ifdef TRANSFER
  if (this->InitializeRadiativeTransferFields() == FAIL) {
    fprintf(stderr, "\nError in InitializeRadiativeTransferFields.\n");
    return FAIL;
  }
#endif

  printf("rho_sphere=%"GSYM", cs_sphere=%"GSYM", rho_medium=%"GSYM", p_medium=%"GSYM"\n",
	 rho_sphere[0], cs_sphere[0], rho_medium, p_medium);

  // if use BE sphere, read in the BE sphere density profile

  char *filename = "be.dat";
  int n_bin = 6401;
  float radius[n_bin];
  float rho_be[n_bin];

  if (sphere_type[0] == 3 || sphere_type[0] == 4) {
    FILE *fptr = fopen(filename, "r");
    char line[MAX_LINE_LENGTH];
    for (int i = 0; i < n_bin; i++) {
      if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL) {
        printf("BE sphere data not enough\n");
        return FAIL;
      }
      sscanf(line, "%"GSYM" %"GSYM, &radius[i], &rho_be[i]);
    }
    fclose(fptr);
  }

  
  float rho, vel[3], eint, etot, h, cs, dpdrho, dpde, v2;
  FLOAT sinphi, cosphi;
  int n = 0;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

	FLOAT x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	FLOAT y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	FLOAT z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	rho = rho_medium;
	EOS(p_medium, rho_medium, eint, h, cs, dpdrho, dpde, EOSType, 1);
	for (dim = 0; dim < 3; dim++) {
	  vel[dim] = 0.0;
	}

	/* Loop over spheres. */
	for (sphere = 0; sphere < n_sphere; sphere++) {
          
	  /* Find distance from center. */

	  FLOAT r = sqrt(pow(fabs(x-sphere_position[sphere][0]), 2) +
		   pow(fabs(y-sphere_position[sphere][1]), 2) +
		   pow(fabs(z-sphere_position[sphere][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);

	  if (r < r_sphere[sphere]) {

            FLOAT xpos, ypos, zpos, drad;
                                                                                                                                                          
	    /* Compute position. */
	    
	    xpos = x-sphere_position[sphere][0];
	    ypos = y-sphere_position[sphere][1];
	    zpos = z-sphere_position[sphere][2];

	    // compute the azimuthal angle
	    cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
	    sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);
	    
	    /* 1. Uniform density */

            if (sphere_type[sphere] == 1) {
              rho = rho_sphere[sphere];
	      EOS(p_sphere[sphere], rho, eint, h, cs, dpdrho, dpde, EOSType, 1);
	      /*float ps;
	      EOS(ps, rho, eint, 0, 2);
	      printf("ps=%"GSYM", pm=%"GSYM", eint=%"GSYM"\n", ps, p_medium, eint);*/
            }
	    
	    /* 2. Uniform, uniformly rotating */

	    if (sphere_type[sphere] == 2) {
	      rho = rho_sphere[sphere];
	      FLOAT cos2phi = cosphi*cosphi -sinphi*sinphi;
	      FLOAT cos3phi = 4.0*pow(cosphi,3) - 3.0*cosphi;
	      //rho *= (1.0 + 0.5*cos2phi);
	      //rho *= (1.0 + 0.5*cos3phi);
	      //EOS(p_sphere[sphere], rho, eint, h, cs, dpdrho, dpde, EOSType, 1);
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	      printf("cs=%"GSYM", vel[0]=%"GSYM" ", cs, vel[0]);
	      double mach_turb = 0.0;
	      vel[0] += mach_turb*Gaussian(cs);
	      vel[1] += mach_turb*Gaussian(cs);
	      vel[2] += mach_turb*Gaussian(cs);
	      //printf("vel[0]=%"GSYM" \n", vel[0]);
	    }

	    /* 3. Bonnor-Ebert sphere */

	    if (sphere_type[sphere] == 3) {
	      double ksi_e = 6.451; // critical radius of BE sphere
	      FLOAT r_be = r*ksi_e/r_sphere[sphere];
	      // find the position of r_be in rho_be
	      FLOAT dis_old = 1e10, dis;
	      int n;
	      for (n = 0; n < n_bin; n++) {
		dis = fabs(radius[n]-r_be);
		if (dis > dis_old) {
		  break;
		} else {
		  dis_old = dis;
		}
	      }
	      if (n == n_bin) {
		n = n_bin -1;
	      }
	      rho = rho_sphere[sphere]*rho_be[n];
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      FLOAT cos2phi = cosphi*cosphi - sinphi*sinphi;
	      FLOAT cos3phi = 4.0*pow(cosphi,3) - 3.0*cosphi;
	      FLOAT sin3phi = 3.0*sinphi - 4.0*pow(sinphi, 3);
	      FLOAT phi0 = M_PI/4.0;
	      FLOAT cos3phi0 = cos3phi*cos(phi0) - sin3phi*sin(phi0);
	      FLOAT omega2 = 0.0;//omega_sphere[sphere]; // velocity perturbation
	      FLOAT omega = omega_sphere[sphere] + omega2*cos2phi;
	      //rho *= (1.0 + 0.5*(cos2phi+cos3phi0));
	      //rho *= (1.0 + 0.1*cos2phi);
	      vel[0] = -omega*ypos;
	      vel[1] = omega*xpos;
	      double mach_turb = 0.0;
	      vel[0] += mach_turb*Gaussian(cs_sphere[sphere]);
	      vel[1] += mach_turb*Gaussian(cs_sphere[sphere]);
	      vel[2] += mach_turb*Gaussian(cs_sphere[sphere]);
	    }	      

	    /* 4. Flattened 1/r^2 sphere */

	    if (sphere_type[sphere] == 4) {
	      rho = rho_sphere[sphere] / (1.0 + pow(3.0*r/r_sphere[sphere],2));
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	    }

	  } // if (r < r_sphere)
	} // end: loop over spheres

	v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
	BaryonField[iden ][n] = rho;
	BaryonField[ivx  ][n] = vel[0];
	BaryonField[ivy  ][n] = vel[1];
	BaryonField[ivz  ][n] = vel[2];
	BaryonField[ietot][n] = eint + 0.5*v2;
	
	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}

        if (MultiSpecies > 0) {
          BaryonField[HIINum][n] = InitialFractionHII *
            CoolData.HydrogenFractionByMass * BaryonField[0][n];
          BaryonField[HeIINum][n] = InitialFractionHeII*
            BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
          BaryonField[HeIIINum][n] = InitialFractionHeIII*
            BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
          BaryonField[HeINum][n] =
            (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
            BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
        
	  BaryonField[HINum][n] =
	    CoolData.HydrogenFractionByMass*BaryonField[0][n] - BaryonField[HIINum][n];

	  if (MultiSpecies > 1) {
	    BaryonField[HINum][n] -= 
	      BaryonField[HMNum][n] + BaryonField[H2IINum][n] + BaryonField[H2INum][n];
	  }

	  // electron "density": n_e * m_p                                                                                     
	  BaryonField[DeNum][n] = BaryonField[HIINum][n] +
	    0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	  if (MultiSpecies > 1) {
	    BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] -
	      BaryonField[HMNum][n];
	  }

	  /* Set Deuterium species (assumed to be negligible). */

	  if (MultiSpecies > 2) {
	    BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
            BaryonField[HINum][n];
	    BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
	      BaryonField[HIINum][n];
	    BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
	      BaryonField[H2INum][n];
	  }
	} // if (MultiSpecies > 0)

      } // end loop over grid
    }
  }

  /* If needed, set a radiation field in the cell where the                                                                  
     sources resides to flag the cells by optical depth */

  /*FLOAT pos[3];
  int ci[3];
  RadiationSourceEntry *Source = GlobalRadiationSources->NextSource;

  if (RefineByOpticalDepth) {

    while (Source != NULL) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
        pos[dim] = Source->Position[dim];
        ci[dim] = int((pos[dim] - GridLeftEdge[dim]) / CellWidth[dim][0]);
      }

      printf("PhotonTest[%"ISYM"]: Left edge = %"FSYM" %"FSYM" %"FSYM"\n", level, GridLeftEdge[0],
             GridLeftEdge[1], GridLeftEdge[2]);
      printf("PhotonTest[%"ISYM"]: source (%"FSYM" %"FSYM" %"FSYM") in %"ISYM" %"ISYM" %"ISYM"\n",
             level, pos[0], pos[1], pos[2], ci[0], ci[1], ci[2]);

      if (pos[0] >= GridLeftEdge[0] && pos[0] <= GridRightEdge[0] &&
          pos[1] >= GridLeftEdge[1] && pos[1] <= GridRightEdge[1] &&
          pos[2] >= GridLeftEdge[2] && pos[2] <= GridRightEdge[2]) {

        index = GRIDINDEX(ci[0], ci[1], ci[2]);
        BaryonField[kphHINum][index] = 1;

        printf("PhotonTest[%"ISYM"]: set kphHI in %"ISYM" %"ISYM" %"ISYM" (%"ISYM")\n",
               level, ci[0], ci[1], ci[2], index);
      } 

      Source = Source->NextSource;

    } 

  } */
  int SetStarParticle = 0;
  
  if (SetStarParticle && level == 0) {
    double mass_p = 10.0 * 1.989e33;
    mass_p /= massu;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*rhou));
    t_dyn /= tu;
      

    NumberOfParticles = 1;
    NumberOfStarParticles = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = 0.501; // 0.6; // 0.55;                                                    
    ParticlePosition[1][0] = 0.501;
    ParticlePosition[2][0] = 0.501;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    ParticleAttribute[0][0] = 0.001; // creation time
    ParticleAttribute[1][0] = t_dyn; // dynamical time
    ParticleAttribute[2][0] = 0.0; // 
    printf("NumberOfParticles = %"ISYM"\n", NumberOfParticles);
  }

  printf("np=%"ISYM"\n", NumberOfParticles);

  return SUCCESS;
}

/* generate random number using gaussian distribution, based on Numerical Recipe */

double Gaussian(double cs)
{

  double mean = 0;
  double stdev = cs;
  double u1 = rand();
  u1 = u1/RAND_MAX;
  double u2 = rand();
  u2 = u2/RAND_MAX;
  double x = mean + stdev*sqrt(-2*log(u1))*cos(2*M_PI*u2);
  
  return x;
}
