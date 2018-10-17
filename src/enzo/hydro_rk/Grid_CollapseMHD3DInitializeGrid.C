/***********************************************************************
/
/  GRID CLASS (INITIALIZE MAGNETIZED CLOUD)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1: Tom Abel 2010, add turbulence generator to this problem as well
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
void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
			  float kmin, float kmax, float dk,
			  FLOAT **LeftEdge, FLOAT **CellWidth, int seed);


int grid::CollapseMHD3DInitializeGrid(int n_sphere,
				      FLOAT r_sphere[MAX_SPHERES],
				      FLOAT rc_sphere[MAX_SPHERES],
				      float rho_sphere[MAX_SPHERES],
				      float p_sphere[MAX_SPHERES],
				      float cs_sphere[MAX_SPHERES],
				      FLOAT sphere_position[MAX_SPHERES][MAX_DIMENSION],
				      float omega_sphere[MAX_SPHERES], 
				      float turb_sphere[MAX_SPHERES], 
				      float Bnaught, float theta_B,
				      int Bdirection,
				      int   sphere_type[MAX_SPHERES],
				      float rho_medium, float p_medium, int level,
				      int SetBaryonFields)
{

  const float HIIFraction = 1.2e-5;
  const float HeIIFraction = 1.0e-14;
  const float HeIIIFraction = 1.0e-17;
  const float HMFraction = 2.0e-9;
  const float H2IFraction = 2.0e-4;
  const float H2IIFraction = 3.0e-14;
  int TurbulenceSeed = 191105;
  float *TurbulenceVelocity[3];
  float HIFraction, HeIFraction, eFraction;
  HIFraction = CoolData.HydrogenFractionByMass - HIIFraction;
  if (MultiSpecies > 1)
    HIFraction -= HMFraction + H2IFraction + H2IIFraction;
  HeIFraction = 1.0 - CoolData.HydrogenFractionByMass - 
    HeIIFraction - HeIIIFraction;
  eFraction = HIIFraction + 0.25*HeIIFraction + 0.5*HeIIIFraction;
  if (MultiSpecies > 1)
    eFraction += 0.5*H2IIFraction - HMFraction;


  /* declarations */

  int dim, i, j, k, m, field, sphere, size, activesize;

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
  if ( UseMHD ) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if( HydroMethod == MHD_RK ){
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

  if(UsePoissonDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
    FieldType[NumberOfBaryonFields++] = DebugField;  
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, 
    H2INum, H2IINum, DINum, DIINum, HDINum;

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

  if (WritePotential) {
    FieldType[NumberOfBaryonFields++] = GravPotential;
    //    FieldType[NumberOfBaryonFields++] = AccelerationField1;
    //    FieldType[NumberOfBaryonFields++] = AccelerationField2;
    //    FieldType[NumberOfBaryonFields++] = AccelerationField3;
  }

  FieldType[NumberOfBaryonFields++] = DebugField;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  //  for (i=0; i< NumberOfBaryonFields; i++)    BaryonField[i] = NULL;


  printf("Grid_CollapseMHD3DInitialize: Setting up grid variables.\n");

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
  activesize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  for (dim = 0; dim < GridRank; dim++) {
    TurbulenceVelocity[dim] = new float[activesize];
    for (int n = 0; n < activesize; n++) {
      TurbulenceVelocity[dim][n] = 0.0;
    }
  }

  this->AllocateGrids(); 

  printf("rho_sphere=%"GSYM", cs_sphere=%"GSYM", rho_medium=%"GSYM", p_medium=%"GSYM"\n",
	 rho_sphere[0], cs_sphere[0], rho_medium, p_medium);
  printf("r_sphere: %"GSYM"\n", r_sphere[0]);
  printf("turb_sphere: %"GSYM"\n", turb_sphere[0]);

  // if use BE sphere, read in the BE sphere density profile

  char *filename = "be.dat";
  int n_bin = 6401;
  float radius[n_bin];
  float rho_be[n_bin];

  if (sphere_type[0] == 3 || sphere_type[0] == 4) {
    printf("Opening ./be.dat\n");
    FILE *fptr = fopen(filename, "r");
    char line[MAX_LINE_LENGTH];
    for (int i = 0; i < n_bin; i++) {
      if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL) {
        printf("BE sphere data not enough\n");
        return FAIL;
      }
      sscanf(line, "%"FSYM" %"FSYM, &radius[i], &rho_be[i]);
    }
    printf("Reading ./be.dat finished.\n");
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
  FLOAT phi, theta, x=0., y=0., z=0.;
  int n = 0;
  
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) {

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1) y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2) z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	rho = rho_medium;
	EOS(p_medium, rho_medium, eint, h, cs, dpdrho, dpde, EOSType, 1);
	for (dim = 0; dim < 3; dim++) {
	  vel[dim] = 0.0;
	}

	Bx = 0.0;
	By = 0.0;
	Bz = 0.0;

	switch (Bdirection) {
	case 0:
	  Bx = Bnaught;
	  break;
	case 1:
	  By = Bnaught;
	  break;
	case 2:
	  Bz = Bnaught;
	  break;
	default:
	  ENZO_FAIL("Bdirection must be 0,1,2");
	}

	/* Loop over spheres. */
	for (sphere = 0; sphere < n_sphere; sphere++) {
          
	  /* Find distance from center. */

	  FLOAT r = sqrt(pow(fabs(x-sphere_position[sphere][0]), 2) +
		   pow(fabs(y-sphere_position[sphere][1]), 2) +
		   pow(fabs(z-sphere_position[sphere][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);

	  if (r < r_sphere[sphere]) {

            FLOAT xpos=0., ypos=0., zpos=0., drad;

	    xpos = x-sphere_position[sphere][0];
	    if (GridRank > 1) ypos = y-sphere_position[sphere][1];
	    if (GridRank > 2) zpos = z-sphere_position[sphere][2];

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
	      float p, cs, h, dpdrho, dpde;
	      p = pow(cs_sphere[sphere], 2)*rho/Gamma;
	      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1); 
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
	      float p, cs, h, dpdrho, dpde;
	      p = rho*pow(cs_sphere[sphere], 2)/Gamma;
	      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1); 
	      //	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
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

	    /* Rotating Gaussian of Truelove et al 1997 */

	    if (sphere_type[sphere] == 7) {
	      float m2mode = 1. + 0.1*cos(2.*phi);
	      rho = rho_sphere[sphere] * exp(-pow(r/r_sphere[sphere]/0.58,2));

	      //	      rho *= m2mode;
	      float p, cs, h, dpdrho, dpde;
	      p = rho*pow(cs_sphere[sphere], 2)/Gamma;
	      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1); 
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	    }

	    /* Rotating Gaussian similar to Truelove et al 1997 no m=2 mode*/

	    if (sphere_type[sphere] == 8) {
	      rho = rho_sphere[sphere] * exp(-pow(r/r_sphere[sphere]/0.58,2));
	      float p, cs, h, dpdrho, dpde;
	      p = rho*pow(cs_sphere[sphere], 2)/Gamma;
	      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1); 
	      vel[0] = -omega_sphere[sphere]*ypos;
	      vel[1] = omega_sphere[sphere]*xpos;
	    }

	    if (sphere_type[sphere] == 9) {
	      rho  = rho_sphere[sphere];
	      // Uniform Burkert Bodenheimer test
	      // BUT keep the pressure constant everywhere
	      // to avoid discontinuities at the sphere boundaries
	      float p, cs, h, dpdrho, dpde;
	      p = pow(cs_sphere[sphere], 2)*rho/Gamma;
	      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1); 
	      // for the B-field we put it along x to slow the 
	      // collapse along the z direction
	      Bx = Bnaught;
	      By = 0;
	      Bz = 0;
              vel[0] = -omega_sphere[sphere]*ypos;
              vel[1] = omega_sphere[sphere]*xpos;
	    }

	    if (sphere_type[sphere] == 10) {
	      rho  = rho_sphere[sphere];
	      // Cooling Sphere from Volker Springel
	      // 
	      // 
	      float p, cs, h, dpdrho, dpde, a2, rho0;
	      a2 = PointSourceGravityCoreRadius*PointSourceGravityCoreRadius;
	      rho = rho_sphere[sphere]*a2/(a2+xpos*xpos);
	      eint = .75 * PointSourceGravityConstant;
	      // for the B-field we put it along x to slow the 
	      // collapse along the z direction
	      Bx = Bnaught;
	      By = 0;
	      Bz = 0;
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
	if (UseMHD) {
	  BaryonField[iBx ][n] = Bx;
	  BaryonField[iBy ][n] = By;
	  BaryonField[iBz ][n] = Bz;
    }
    if( HydroMethod == MHD_RK ){
	  BaryonField[iPhi][n] = 0.0;
	}
    if( UseMHDCT ){
        MagneticField[0][n] = Bx;
        MagneticField[1][n] = By;
        MagneticField[2][n] = Bz;
    }
	BaryonField[NumberOfBaryonFields-1][n] = 0.;
      } // end loop over grid
    }
  }



  /* Initialize turbulent velocity field */

  if (SetBaryonFields && (turb_sphere[0] > 0.)) {

    float k1, k2, dk;
      k1 = 5;
      k2 = 8.0;
      dk = 1.0;

    printf("Begin generating turbulent velocity spectrum...\n");
    Turbulence_Generator(TurbulenceVelocity, 
			 GridDimension[0]-2*NumberOfGhostZones, 
			 GridDimension[1]-2*NumberOfGhostZones,
			 GridDimension[2]-2*NumberOfGhostZones,
			 4.0, k1, k2, dk,
			 CellLeftEdge, CellWidth, TurbulenceSeed);    

    printf("Turbulent spectrum generated\n");

    float VelocityNormalization = 1;
// for level > 0 grids the CloudMachNumber passed in is actuall the Velocity normalization factor
    if (level > 0) VelocityNormalization = turb_sphere[0];
    printf("Cloud Mach Number = %"GSYM" \n",turb_sphere[0]);
    for (i = 0; i < 3; i++) {
      for (n = 0; n < activesize; n++) {
	TurbulenceVelocity[i][n] *= VelocityNormalization;
      }
    }


    /* Set turbulent velocity field */
  FLOAT x,y,z;
  int igrid;
    n = 0;
    /*    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) { */
	  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { 
	    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) { 
	  igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  
	  BaryonField[ivx][igrid] += TurbulenceVelocity[0][n];
	  BaryonField[ivy][igrid] += TurbulenceVelocity[1][n];
	  BaryonField[ivz][igrid] += TurbulenceVelocity[2][n];
	  BaryonField[ietot][igrid] += 
	    0.5 * (pow(TurbulenceVelocity[0][n],2) + 
		   pow(TurbulenceVelocity[1][n],2) + 
		   pow(TurbulenceVelocity[2][n],2));

	} 
      }
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
  float ppos[3];
  ppos[0] = 0.50001;
  ppos[1] = 0.50001;
  ppos[2] = 0.50001;

  this->DeleteParticles();
  NumberOfParticles = 0;
  NumberOfStars    = 0;
  //  if (PutSinkParticle == 1 && level == MaximumRefinementLevel) {
  if (PutSinkParticle == 1 && GridLeftEdge[0] < ppos[0] && GridRightEdge[0] > ppos[0] &&
      GridLeftEdge[1] < ppos[1] && GridRightEdge[1] > ppos[1] &&
      GridLeftEdge[2] < ppos[2] && GridRightEdge[2] > ppos[2]) {

    double mass_p = 1.*1.989e32;
    mass_p /= MassUnits;
    double dx = CellWidth[0][0];
    double den_p = mass_p / pow(dx,3);
    double t_dyn = sqrt(3*M_PI/(6.672e-8*den_p*DensityUnits));
    t_dyn /= TimeUnits;


    NumberOfParticles = 1;
    NumberOfStars = 1;
    NumberOfParticleAttributes = 3;
    //    MaximumParticleNumber = 1;
    this->AllocateNewParticles(NumberOfParticles);
    ParticleMass[0] = den_p;
    ParticleNumber[0] = 0;
    ParticleType[0] = PARTICLE_TYPE_MUST_REFINE;
    ParticlePosition[0][0] = ppos[0]; // 0.6; // 0.55;                         
    ParticlePosition[1][0] = ppos[1];
    ParticlePosition[2][0] = ppos[2];

    ParticleVelocity[0][0] = 0.0;
    ParticleVelocity[1][0] = 0.0;
    ParticleVelocity[2][0] = 0.0;
    ParticleAttribute[0][0] = 0.01; // creation time
    ParticleAttribute[1][0] = 0; // dynamical time
    ParticleAttribute[2][0] = mass_p;
  }

  /* Set uniform species fractions */

  int index;
  if (MultiSpecies > 0) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  BaryonField[DeNum][index] = eFraction * BaryonField[0][index];
	  BaryonField[HINum][index] = HIFraction * BaryonField[0][index];
	  BaryonField[HIINum][index] = HIIFraction * BaryonField[0][index];
	  BaryonField[HeINum][index] = HeIFraction * BaryonField[0][index];
	  BaryonField[HeIINum][index] = HeIIFraction * BaryonField[0][index];
	  BaryonField[HeIIINum][index] = HeIIIFraction * BaryonField[0][index];
	}
  }

  if (MultiSpecies > 1) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  BaryonField[HMNum][index] = HMFraction * BaryonField[0][index];
	  BaryonField[H2INum][index] = H2IFraction * BaryonField[0][index];
	  BaryonField[H2IINum][index] = H2IIFraction * BaryonField[0][index];
	}
  }

  if (MultiSpecies > 2) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  BaryonField[DINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[HINum][index];
	  BaryonField[DIINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[HIINum][index];
	  BaryonField[HDINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[H2INum][index];
	}
  }

  printf("Grid_CollapseMHD3DInitiailize: done with this grid\n");
  return SUCCESS;
}
