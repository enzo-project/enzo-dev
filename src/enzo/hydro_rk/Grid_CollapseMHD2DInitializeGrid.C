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

int grid::CollapseMHD2DInitializeGrid(FLOAT r_sphere,
				      FLOAT rc_sphere,
				      float rho_sphere,
				      float p_sphere, 
				      float cs_sphere,
				      float omega_sphere, float B0, float theta_B,
				      int   sphere_type,
				      float rho_medium, float p_medium, int level)
{
  /* declarations */

  int dim, i, j, k, m, sphere;

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

  if(UsePoissonDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
    FieldType[NumberOfBaryonFields++] = DebugField;  
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  /* Units and parameters */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  float MagneticUnits = sqrt(4.0*Pi*DensityUnits)*VelocityUnits;
  double G = 6.67e-8;

  this->AllocateGrids();

  printf("rho_sphere=%"GSYM", cs_sphere=%"GSYM", rho_medium=%"GSYM", p_medium=%"GSYM"\n",
	 rho_sphere[0], cs_sphere[0], rho_medium*DensityUnits, p_medium);

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
	Bz = B0;

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

	    
	    /* 0. uniform sphere */
	    
	    if (sphere_type[sphere] == 0) {
	      rho  = rho_sphere[sphere];
	      eint = pow(cs_sphere[sphere], 2)/(Gamma-1.0);
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

  return SUCCESS;
}
