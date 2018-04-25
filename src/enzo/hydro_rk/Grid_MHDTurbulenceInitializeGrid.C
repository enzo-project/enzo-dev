/***********************************************************************
/
/  GRID CLASS (INITIALIZE MAGNETIZED TURBULENT CLOUD)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1: Tom Abel adopted  (2009)
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

int grid::MHDTurbulenceInitializeGrid(float rho_medium, float cs_medium, float mach, 
				      float Bnaught, int seed, int level, int SetBaryonFields)
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

  int idrivex, idrivey, idrivez;
  if (UseDrivingField && (HydroMethod == HD_RK || HydroMethod == MHD_RK)) {
    idrivex = NumberOfBaryonFields;
    idrivey = idrivex + 1;
    idrivez = idrivex + 2;
    FieldType[NumberOfBaryonFields++] = DrivingField1;
    FieldType[NumberOfBaryonFields++] = DrivingField2;
    FieldType[NumberOfBaryonFields++] = DrivingField3;
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

  this->AllocateGrids();
  
  for (int dim=0; dim < GridRank; dim++) 
    if (UseDrivingField && (HydroMethod < 3) && RandomForcingField[dim] == NULL) {
      fprintf(stderr,"Driving turbulence only implemented for Hydromethod > 2 !\n");
      fprintf(stderr,"Or turn UseDrivingField off to study decaying turbulence \n");
      fprintf(stderr,"with this hydro solver. \n");
      ERROR_MESSAGE;
      //      RandomForcingField[dim] = new float[size];
    }


  int activesize = 1,n=0;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  float *TurbulenceVelocity[3],*DrivingField[3];
  for (int dim = 0; dim < 3; dim++) {
    TurbulenceVelocity[dim] = new float[activesize];
    if (UseDrivingField) DrivingField[dim] = new float[activesize];
    for (n = 0; n < activesize; n++) {
      TurbulenceVelocity[dim][n] = 0.0;
    if (UseDrivingField) DrivingField[dim][n] = 0.0;
    }
  }

  if (debug) 
    printf("Begin generating turbulent velocity spectrum... %"ISYM" %"ISYM" %"ISYM"\n", 
	   GridDimension[0]-2*NumberOfGhostZones,
	   GridDimension[1]-2*NumberOfGhostZones,
	   GridDimension[2]-2*NumberOfGhostZones);

  Turbulence_Generator(TurbulenceVelocity, GridDimension[0]-2*NumberOfGhostZones,
		       GridDimension[1]-2*NumberOfGhostZones,
		       GridDimension[2]-2*NumberOfGhostZones,
		       4.0, 1, 5, 1,
		       CellLeftEdge, CellWidth, seed);
  printf("Turbulent spectrum generated\n");

  float VelocityNormalization = 1;
// for level > 0 grids the CloudMachNumber passed in is actuall the Velocity normalization factor
  if (level > 0) VelocityNormalization = mach; 

  for (int i = 0; i < 3; i++) {
    for (n = 0; n < activesize; n++) {
      TurbulenceVelocity[i][n] *= VelocityNormalization;
    }
  }


  // assume isothermal initially
  float p_medium = rho_medium*cs_medium*cs_medium;


  // For coding efficiency, we set the radius, density contrast of the 
  // initial cloud by hand below. Should be changed in the future
  // Initialize field without turbulent velocity field
  float eint, h, dpdrho, dpde, cs;
  eint = cs_medium*cs_medium/(Gamma-1.0);
  FLOAT xc = 0.5, yc = 0.5, zc = 0.5, x, y, z, r;
  FLOAT rs = 1.; // 0.3;
  n=0;
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
	  BaryonField[ietot][n] = eint + 0.5*Bnaught*Bnaught/rho_medium;
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
	  BaryonField[ietot][n] = eint_out + 0.5*Bnaught*Bnaught/rho_out;
          if (DualEnergyFormalism) {
            BaryonField[ieint][n] = eint_out;
          }
	}

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][n] = 0.0;
	  BaryonField[iBy ][n] = 0.0;
	  BaryonField[iBz ][n] = Bnaught;
	  BaryonField[iPhi][n] = 0.0;
	  //	  BaryonField[ietot][n] += 0.5 * Bnaught*Bnaught / BaryonField[iden ][n] ;
	}
	if (UseDrivingField && (HydroMethod == HD_RK || HydroMethod == MHD_RK)) {
	  BaryonField[idrivex][n] = 0.0;
	  BaryonField[idrivey][n] = 0.0;
	  BaryonField[idrivez][n] = 0.0;
	}
	if (UseDrivingField && (HydroMethod < 3)) {
	  RandomForcingField[0][n] = 0.;
	  RandomForcingField[1][n] = 0.;
	  RandomForcingField[2][n] = 0.;
	}
      }
    }
  }

  // Set the turbulent velocity field
  float v2;

  int igrid;
  n = 0;
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
          BaryonField[ivx][igrid] = TurbulenceVelocity[0][n];
          BaryonField[ivy][igrid] = TurbulenceVelocity[1][n];
          BaryonField[ivz][igrid] = TurbulenceVelocity[2][n];
          BaryonField[ietot][igrid] += 
	    0.5*(TurbulenceVelocity[0][n]*TurbulenceVelocity[0][n] + 
		 TurbulenceVelocity[1][n]*TurbulenceVelocity[1][n] + 
		 TurbulenceVelocity[2][n]*TurbulenceVelocity[2][n]);
        }

      } 
    }
  }

  for (int i = 0; i < 3; i++) {
    delete [] TurbulenceVelocity[i];
    }

  /* Initialize driving force field = efficiency * density * velocity / t_ff*/
  printf("UseDrivingField =%"ISYM"\n",UseDrivingField);
  if (UseDrivingField) {
    float k1, k2, dk;
    k1 = 3.0;
    k2 = 4.0;
    dk = 1.0;
    printf("Begin generating driving force field ...\n");
    Turbulence_Generator(DrivingField, GridDimension[0]-2*NumberOfGhostZones, 
		       GridDimension[1]-2*NumberOfGhostZones,
		       GridDimension[2]-2*NumberOfGhostZones, 
			 4.0, k1, k2, dk,
			 CellLeftEdge, CellWidth, seed);
    printf("Driving force field generated\n");


    /* Renormalize to ensure <F>=0 */

    double Fx = 0.0, Fy = 0.0, Fz = 0.0;
    for (n = 0; n < activesize; n++) {
      Fx += DrivingField[0][n];
      Fy += DrivingField[1][n];
      Fz += DrivingField[2][n];
    }

    Fx /= activesize;
    Fy /= activesize;
    Fz /= activesize;

    for (n = 0; n < activesize; n++) {
      DrivingField[0][n] -= Fx;
      DrivingField[1][n] -= Fy;
      DrivingField[2][n] -= Fz;
    }


    n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	  igrid = i + GridDimension[0]*(j+k*GridDimension[1]);
	  if (HydroMethod == HD_RK || HydroMethod == MHD_RK) {
	    BaryonField[idrivex][igrid] = DrivingField[0][n];
	    BaryonField[idrivey][igrid] = DrivingField[1][n];
	    BaryonField[idrivez][igrid] = DrivingField[2][n];
	  } else {
	     RandomForcingField[0][igrid] = DrivingField[0][n]*DrivingEfficiency;
	     RandomForcingField[1][igrid] = DrivingField[1][n]*DrivingEfficiency;
	     RandomForcingField[2][igrid] = DrivingField[2][n]*DrivingEfficiency;
	     //	     fprintf(stderr, "%"GSYM"\t",RandomForcingField[0][igrid]);
	  }
	}
      }
    }


     for (int dim = 0; dim < GridRank; dim++) {
      delete [] DrivingField[dim];
      }
  }    

  //printf("Grid_MHDTurb: COMPLETED\n");
  return SUCCESS;
}

