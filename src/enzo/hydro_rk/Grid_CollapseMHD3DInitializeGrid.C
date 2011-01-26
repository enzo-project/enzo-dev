/***********************************************************************
/
/  GRID CLASS (INITIALIZE MAGNETIZED CLOUD)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
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
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double Gaussian(double cs);

int grid::CollapseMHD3DInitializeGrid(int n_sphere,
				      FLOAT r_sphere[MAX_SPHERES],
				      FLOAT rc_sphere[MAX_SPHERES],
				      float rho_sphere[MAX_SPHERES],
				      float p_sphere[MAX_SPHERES],
				      float cs_sphere[MAX_SPHERES],
				      FLOAT sphere_position[MAX_SPHERES][MAX_DIMENSION],
				      float omega_sphere[MAX_SPHERES], float Bnaught, float theta_B,
				      int   sphere_type[MAX_SPHERES],
				      float rho_medium, float p_medium, int level)
{
  /* declarations */

  int dim, i, j, k, m, field, sphere, size;

    int phip_num;
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

  if(UseDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
    FieldType[NumberOfBaryonFields++] = DebugField;  
  }

  if (WritePotential) {
    FieldType[NumberOfBaryonFields++] = GravPotential;
    FieldType[NumberOfBaryonFields++] = AccelerationField1;
    FieldType[NumberOfBaryonFields++] = AccelerationField2;
    FieldType[NumberOfBaryonFields++] = AccelerationField3;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  /* Units and parameters */

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  float MagneticUnits = sqrt(4.0*M_PI*DensityUnits)*VelocityUnits;
  double G = 6.67e-8;

  size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int count=0;
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
      count++;
    }
  }
  printf("Allocated %"ISYM" Baryonfields\n", count);

  printf("rho_sphere=%"GSYM", cs_sphere=%"GSYM", rho_medium=%"GSYM", p_medium=%"GSYM"\n",
	 rho_sphere[0], cs_sphere[0], rho_medium, p_medium);
  printf("r_sphere: %"GSYM"\n", r_sphere[0]);
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

  /* If use SIT, read in the SIT sphere density profile */

  float theta_sit[n_bin];
  float R_sit[n_bin], phi_sit[n_bin], dphi_sit[n_bin], v_sit;

  if (sphere_type[0] == 5) {
    filename = "sit.dat";
    n_bin = 1000;
    FILE *fptr = fopen(filename, "r");
    char line[MAX_LINE_LENGTH];
    /* Get rid of the comment line */
    fgets(line, MAX_LINE_LENGTH, fptr);
    for (int i = 0; i < n_bin; i++) {
      if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL) {
        printf("SIT data not enough\n");
        return FAIL;
      }
      sscanf(line, "%"GSYM" %"GSYM" %"GSYM" %"GSYM, &theta_sit[i], &R_sit[i], &phi_sit[i], &dphi_sit[i]);
    }
    fgets(line, MAX_LINE_LENGTH, fptr);
    sscanf(line, "%"GSYM, &v_sit);
    fclose(fptr);
  }

  float rho, vel[3], eint, etot, h, cs, dpdrho, dpde, v2, B2, Bx, By, Bz;
  FLOAT phi, theta;
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
	Bx = 0.0;
	By = 0.0;
	Bz = Bnaught;

	/* Loop over spheres. */
	for (sphere = 0; sphere < n_sphere; sphere++) {
          
	  /* Find distance from center. */

	  FLOAT r = sqrt(pow(fabs(x-sphere_position[sphere][0]), 2) +
		   pow(fabs(y-sphere_position[sphere][1]), 2) +
		   pow(fabs(z-sphere_position[sphere][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);

	  if (r < r_sphere[sphere]) {

            FLOAT xpos, ypos, zpos, drad;

	    xpos = x-sphere_position[sphere][0];
	    ypos = y-sphere_position[sphere][1];
	    zpos = z-sphere_position[sphere][2];

	    // compute the azimuthal angle
	    FLOAT cosphi = xpos/sqrt(xpos*xpos+ypos*ypos);
	    FLOAT sinphi = ypos/sqrt(xpos*xpos+ypos*ypos);

	    /* Compute the azimuthal and polar angles */
	    phi   = acos(xpos/sqrt(xpos*xpos+ypos*ypos));
	    if (ypos < 0) phi = 2.0*M_PI-phi;
	    FLOAT R1 = sqrt(xpos*xpos+ypos*ypos);
	    theta = acos(zpos/r);
	    /*if (fabs(zpos) < 1e-3) {
	      theta = acos(zpos/r);
	    } else if (R1 > CellWidth[0][0]) {
	      theta = acos(zpos/r);
	    } else {
	      theta = acos(zpos/sqrt(pow(xpos+sign(xpos)*CellWidth[0][0],2)+ypos*ypos+zpos*zpos));
	      }*/
	    
	    /* 0. uniform sphere */
	    
	    if (sphere_type[sphere] == 0) { // set field along x for this problem
	      Bx = Bnaught;
	      By = 0.0;
	      Bz = 0.0;
	    };


	    if (sphere_type[sphere] == 0) {
	      rho  = rho_sphere[sphere];
	      FLOAT cos2phi = cosphi*cosphi -sinphi*sinphi;
	      //	      rho *= (1.0 + 0.2*cos2phi);
	      // Burkert & Bodenheimer (1993) m=2 perturbation: 	      
	      float m2mode = 1. + 0.1*cos(2.*phi);
	      rho *= m2mode;
	      // BUT keep the pressure constant everywhere
	      // to avoid discontinuities at the sphere boundaries
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0)/m2mode;
	      // for the B-field we put it along x to slow the 
	      // collapse along the z direction
	      Bx = Bnaught;
	      By = 0;
	      Bz = 0;
              vel[0] = -omega_sphere[sphere]*ypos;
              vel[1] = omega_sphere[sphere]*xpos;
	    }

	    /* 1. Flattened 1/r^2 sphere */

	    if (sphere_type[sphere] == 1) {
	      rho = rho_sphere[sphere] / (1.0 + pow(3.0*r/r_sphere[sphere],2));
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	    }

	    /* 2. Singular Isothermal Sphere */

	    if (sphere_type[sphere] == 2) {
	      rho = pow(cs_sphere[sphere]*VelocityUnits,2)/(2.0*M_PI*G*pow(r*LengthUnits,2));
	      rho /= DensityUnits;
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	    }

	    /* 3. Bonner-Ebert sphere */

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
	      FLOAT omega = omega_sphere[sphere];
	      vel[0] = -omega*ypos;
	      vel[1] = omega*xpos;
	      double mach_turb = 0.0;
	      vel[0] += mach_turb*Gaussian(cs_sphere[sphere]);
	      vel[1] += mach_turb*Gaussian(cs_sphere[sphere]);
	      vel[2] += mach_turb*Gaussian(cs_sphere[sphere]);
	    }	      

	    /* 5. Singular Isothermal Toroid */

	    if (sphere_type[sphere] == 5) {
	      float Br, Btheta;
	      FLOAT dis_theta_old = 1e10, dis_theta;
	      int ii;
	      float theta1 = (theta < M_PI/2.0) ? theta : M_PI - theta;
	      for (ii = 0; ii < n_bin; ii++) {
		dis_theta = fabs(theta_sit[ii] - theta1);
		if (dis_theta > dis_theta_old) {
		  break;
		} else {
		  dis_theta_old = dis_theta;
		}
	      }
	      if (ii == n_bin) {
		ii = n_bin -1;
	      }
	      rho = pow(cs_sphere[sphere]*VelocityUnits,2)*R_sit[ii] / 
		(2.0*M_PI*G*pow(r*LengthUnits,2)) / DensityUnits;
	      Br = 2.0*pow(cs_sphere[sphere]*VelocityUnits,2)*dphi_sit[ii] /
		(r*LengthUnits*sin(theta)*pow(G,0.5)) / MagneticUnits;
	      Btheta = -2.0*pow(cs_sphere[sphere]*VelocityUnits,2)*phi_sit[ii] /
		(r*LengthUnits*sin(theta)*pow(G,0.5)) / MagneticUnits;
	      if (theta > M_PI/2.0) Br *= -1.0;
	      Bx = Br*sin(theta)*cos(phi) + Btheta*cos(theta)*cos(phi);
	      By = Br*sin(theta)*sin(phi) + Btheta*cos(theta)*sin(phi);
	      Bz = Br*cos(theta) - Btheta*sin(theta);
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      vel[0] = -cs_sphere[sphere]*v_sit*sin(phi);
	      vel[1] =  cs_sphere[sphere]*v_sit*cos(phi);
	    }	      

	    /* 1. Flattened 1/r^2 sphere */

	    if (sphere_type[sphere] == 6) {
	      rho = rho_sphere[sphere] / (1.0 + pow(3.0*r/r_sphere[sphere],2));
	      if (theta > M_PI/2.0) theta = M_PI-theta;
	      rho *= (9.0/10.0*theta/(M_PI/2.0) + 1.0/10.0);
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	    }


	  } // if (r < r_sphere)
	} // end: loop over spheres

	v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
	B2 = Bx*Bx + By*By + Bz*Bz;
	BaryonField[iden ][n] = rho;
	BaryonField[ivx  ][n] = vel[0];
	BaryonField[ivy  ][n] = vel[1];
	BaryonField[ivz  ][n] = vel[2];
	BaryonField[ietot][n] = eint + 0.5*v2 + 0.5*B2/rho;
	
	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][n] = Bx;
	  BaryonField[iBy ][n] = By;
	  BaryonField[iBz ][n] = Bz;
	  BaryonField[iPhi][n] = 0.0;
	}

      } // end loop over grid
    }
  }

  int TestBinary = 0;
  if (TestBinary == 1 && level == 0) {

    //NumberOfParticleAttributes = 3;

    double mass_p = 1.0;
    double a = 0.2;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3 * M_PI / den_p);


    NumberOfParticles = 2;
    //NumberOfStarParticles = 2;
    //MaximumParticleNumber = 2;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_DARK_MATTER;
    ParticlePosition[0][0] = 0.5 + 0.5*a;
    ParticlePosition[1][0] = 0.5;
    ParticlePosition[2][0] = 0.5;

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = -sqrt(mass_p/(2.0*a));
    ParticleVelocity[2][0] = 0.0;

    ParticleMass[1] = den_p;
    ParticleNumber[1] = 1;
    ParticleType[1] = PARTICLE_TYPE_DARK_MATTER;
    ParticlePosition[0][1] = 0.5 - 0.5*a;
    ParticlePosition[1][1] = 0.5;
    ParticlePosition[2][1] = 0.5;

    ParticleVelocity[0][1] = 0.0;
    ParticleVelocity[1][1] = sqrt(mass_p/(2.0*a));
    ParticleVelocity[2][1] = 0.0;

  }


  int PutSinkParticle = 0;
  if (PutSinkParticle == 1 && level == MaximumRefinementLevel) {

    double mass_p = 1.1*1.989e33;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;


    NumberOfParticles = 1;
    NumberOfStars = 1;
    //    MaximumParticleNumber = 1;
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
    ParticleAttribute[0][0] = 0.01; // creation time    
    ParticleAttribute[1][0] = 0; // dynamical time                                                                
    ParticleAttribute[2][0] = mass_p;
  }


  return SUCCESS;
}
