/***********************************************************************
/
/  GRID CLASS (INITIALIZE MAGNETIZED TURBULENT CLOUD)
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

void Turbulence_Generator(float **vel, int size, int ind, float sigma, float kmin, float kmax, float dk,
			  FLOAT **LeftEdge, FLOAT **CellWidth, int seed, int level);

int grid::MHDTurbulenceInitializeGrid(float rho_medium, float cs_medium, float mach, 
				      float B0, int seed, int level)
{

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  FieldType[NumberOfBaryonFields++] = Bfield1;
  FieldType[NumberOfBaryonFields++] = Bfield2;
  FieldType[NumberOfBaryonFields++] = Bfield3;
  FieldType[NumberOfBaryonFields++] = PhiField;

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float rhou, lenu, tempu, tu,
    velu, CriticalDensity = 1, BoxLength = 1;
  
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
  
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
    }
  }
  
  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);
  }

  float *vel[3];
  for (int i = 0; i < 3; i++) {
    vel[i] = new float[activesize];
  }

  printf("Begin generating turbulent velocity spectrum...\n");
  Turbulence_Generator(vel, GridDimension[0]-2*DEFAULT_GHOST_ZONES, 4.0, cs_medium*mach, 1, 5, 1,
		       CellLeftEdge, CellWidth, seed, level);
  printf("Turbulent spectrum generated\n");

  // assume isothermal initially
  float p_medium = rho_medium*cs_medium*cs_medium;


  // For coding efficiency, we set the radius, density contrast of the 
  // initial cloud by hand below. Should be changed in the future
  // Initialize field without turbulent velocity field
  float eint, h, dpdrho, dpde, cs;
  eint = cs_medium*cs_medium/(Gamma-1.0);
  int n = 0;
  FLOAT xc = 0.5, yc = 0.5, zc = 0.5, x, y, z, r;
  FLOAT rs = 2.0;
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++, n++) {

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
        r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
        r = max(r, 0.1*CellWidth[0][0]);

        if (r < rs) {
          BaryonField[iden ][n] = rho_medium;
          BaryonField[ivx  ][n] = 0.0;
          BaryonField[ivy  ][n] = 0.0;
          BaryonField[ivz  ][n] = 0.0;
          BaryonField[ietot][n] = eint + 0.5*B0*B0/rho_medium;
          if (DualEnergyFormalism) {
            BaryonField[ieint][n] = eint;
          }
        }
	else {
	  float rho_out = rho_medium/10.0;
	  float eint_out = eint*10.0;
	  BaryonField[iden ][n] = rho_out;
          BaryonField[ivx  ][n] = 0.0;
          BaryonField[ivy  ][n] = 0.0;
          BaryonField[ivz  ][n] = 0.0;
          BaryonField[ietot][n] = eint_out + 0.5*B0*B0/rho_out;
          if (DualEnergyFormalism) {
            BaryonField[ieint][n] = eint_out;
          }
	}
	BaryonField[iBx ][n] = 0.0;
	BaryonField[iBy ][n] = 0.0;
	BaryonField[iBz ][n] = B0;
	BaryonField[iPhi][n] = 0.0;
      }
    }
  }

  // Set the turbulent velocity field
  float v2;
  n = 0;
  int igrid;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {

	igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
        r = sqrt(pow(fabs(x-xc),2)+pow(fabs(y-yc),2)+pow(fabs(z-zc),2));
        r = max(r, 0.1*CellWidth[0][0]);

        if (r < rs) {
          BaryonField[ivx][igrid] = vel[0][n];
          BaryonField[ivy][igrid] = vel[1][n];
          BaryonField[ivz][igrid] = vel[2][n];
          BaryonField[ietot][igrid] += 
	    0.5*(vel[0][n]*vel[0][n] + vel[1][n]*vel[1][n] + vel[2][n]*vel[2][n]);
        }

      } 
    }
  }

  for (int i = 0; i < 3; i++) {
    delete vel[i];
  }

  return SUCCESS;
}

