/***********************************************************************
/
/  GRID CLASS (INITIALIZE RANDOM MAGNETIC FIELD )
/
/  written by: Tom Abel
/  date:       June, 2011, adopted form MHDTurbulenceInitializeGrid
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

void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind, 
			  float kmin, float kmax, float dk,
			  FLOAT **LeftEdge, FLOAT **CellWidth, int seed);

int grid::MHDDecayingRandomFieldInitializeGrid(float rho_medium, float cs_medium, float mach, 
					       float Bnaught, int seed, 
					       float Sindex, float Skmin, float Skmax, 
					       int level, int SetBaryonFields)
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
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }


  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;
  
  if (SetBaryonFields == 0) 
    return SUCCESS;
  

  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0,
    velu = 1.0, CriticalDensity = 1, BoxLength = 1;
  if (UsePhysicalUnit)
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
  

  int activesize = 1,n=0;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  float *RandomBField[3];
  for (int dim = 0; dim < 3; dim++) {
    RandomBField[dim] = new float[activesize];
    for (n = 0; n < activesize; n++) {
      RandomBField[dim][n] = 0.0;
    }
  }

  if (debug && (MyProcessorNumber == ROOT_PROCESSOR)) 
    printf("Begin generating random magnetic field  spectrum... %"ISYM" %"ISYM" %"ISYM"\n", 
	   GridDimension[0]-2*NumberOfGhostZones,
	   GridDimension[1]-2*NumberOfGhostZones,
	   GridDimension[2]-2*NumberOfGhostZones);

  Turbulence_Generator(RandomBField, GridDimension[0]-2*NumberOfGhostZones,
		       GridDimension[1]-2*NumberOfGhostZones,
		       GridDimension[2]-2*NumberOfGhostZones,
		       Sindex, Skmin, Skmax, 1,
		       CellLeftEdge, CellWidth, seed);
  if (MyProcessorNumber == ROOT_PROCESSOR) printf("Random B field spectrum generated\n");

  float BFieldNormalization = 1;
// for level > 0 grids the MachNumber passed in is actuall the Bfield normalization factor
  if (level > 0) BFieldNormalization = mach; 

  for (int i = 0; i < 3; i++) {
    for (n = 0; n < activesize; n++) {
      RandomBField[i][n] *= BFieldNormalization;
    }
  }


  // assume isothermal initially
  float p_medium = rho_medium*cs_medium*cs_medium;


  // For coding efficiency, we set the radius, density contrast of the 
  // initial cloud by hand below. Should be changed in the future
  // Initialize field without turbulent velocity field
  float eint, h, dpdrho, dpde, cs;
  eint = cs_medium*cs_medium;
  n=0;
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++, n++) {

	BaryonField[iden ][n] = rho_medium;
	BaryonField[ivx  ][n] = 0.0;
	BaryonField[ivy  ][n] = 0.0;
	BaryonField[ivz  ][n] = 0.0;
	BaryonField[ietot][n] = eint; // + 0.5*(BaryonField[ivx  ][n]*BaryonField[ivx  ][n]);
	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}
	
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][n] = 0.0;
	  BaryonField[iBy ][n] = 0.0;
	  BaryonField[iBz ][n] = Bnaught;
	  BaryonField[iPhi][n] = 0.0;
	}
      }
    }
  }

  // Set the Random magnetic field
  int igrid;
  n = 0;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	igrid = i + GridDimension[0]*(j+k*GridDimension[1]);

          BaryonField[iBx][igrid] += RandomBField[0][n];
          BaryonField[iBy][igrid] += RandomBField[1][n];
          BaryonField[iBz][igrid] += RandomBField[2][n];
          BaryonField[ietot][igrid] += 
	    0.5*(RandomBField[0][n]*RandomBField[0][n] + 
		 RandomBField[1][n]*RandomBField[1][n] + 
		 RandomBField[2][n]*RandomBField[2][n])/BaryonField[iden ][igrid];

      } 
    }
  }

  for (int i = 0; i < 3; i++) {
    delete [] RandomBField[i];
    }


  //printf("Grid_MHDDecayingRandomField: COMPLETED\n");
  return SUCCESS;
}

