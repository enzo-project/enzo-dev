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

int grid::Collapse1DInitializeGrid(FLOAT r_sphere,
				   FLOAT rc_sphere,
				   float rho_sphere,
				   float p_sphere,
				   float cs_sphere,
				   float omega_sphere,
				   int   sphere_type,
				   float rho_medium, float p_medium)
{
  /* declarations */

  int dim, i, j, k, m, sphere;

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }


  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }


  float rhou, lenu, tempu, tu, velu;
  
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);

  this->AllocateGrids();

  // if use BE sphere, read in the BE sphere density profile

  char *filename = "be.dat";
  int n_bin = 6401;
  float radius[n_bin];
  float rho_be[n_bin];

  if (sphere_type == 3) {
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


  printf("rho_sphere=%"GSYM", cs_sphere=%"GSYM", rho_medium=%"GSYM", p_medium=%"GSYM"\n",
	 rho_sphere, cs_sphere, rho_medium, p_medium);
  
  float rho, vel[3], eint, etot, h, cs, dpdrho, dpde, v2;
  for (i = 0; i < GridDimension[0]; i++) {
    
    FLOAT r = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

    rho = rho_medium;
    EOS(p_medium, rho_medium, eint, h, cs, dpdrho, dpde, 0, 1);
    for (dim = 0; dim < 3; dim++) {
	  vel[dim] = 0.0;
    }

    if (r < r_sphere) {

      if (sphere_type == 1) { // 1. Uniform density
	rho = rho_sphere;
	EOS(p_sphere, rho, eint, h, cs, dpdrho, dpde, 0, 1);
      }
	    
      if (sphere_type == 2) { // 2. Uniform, uniformly rotating
	rho = rho_sphere;
	EOS(p_sphere, rho, eint, h, cs, dpdrho, dpde, 0, 1);
	vel[2] = omega_sphere*r;
      }

      if (sphere_type == 3) { // 3. Bonnor-Ebert sphere
	double ksi_e = 6.451; // critical radius of BE sphere
	FLOAT r_be = r*ksi_e/r_sphere;
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
	rho = rho_sphere*rho_be[n];
	eint = pow(cs_sphere, 2)/(Gamma-1.0);
      }	      
    } // if (r < r_sphere)
    /*else {
      if (sphere_type == 3) {
	rho = max(rho_sphere/14.0*exp(-10.0*(r-r_sphere)/r_sphere), rho_medium);
	EOS(p_medium, rho, eint, h, cs, dpdrho, dpde, 0, 1);
      }
      }*/

    v2 = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
    BaryonField[iden ][i] = rho;
    BaryonField[ivx  ][i] = vel[0];
    BaryonField[ivy  ][i] = vel[1];
    BaryonField[ivz  ][i] = vel[2];
    BaryonField[ietot][i] = eint + 0.5*v2;
    if (DualEnergyFormalism) {
      BaryonField[ieint][i] = eint;
    }
    
  } // end loop over grid

  return SUCCESS;
}

