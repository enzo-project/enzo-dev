/***********************************************************************
/
/  GRID CLASS (INITIALIZE SHEARING BOX)
/
/  written by: Fen Zhao
/  date:       2008
/  modified1: Peng Wang
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
double Gaussian(double cs);

int grid::ShearingBoxInitializeGrid(float AngularVelocity, float VelocityGradient, float ThermalMagneticRatio, 
				    int ShearingBoxProblemType)
{
  /* declarations */


  int phip_num;
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[ieint=NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  
  if(UseDivergenceCleaning) {
    FieldType[NumberOfBaryonFields++] = Phi_pField;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
      for (int i = 0; i < size; i++)
	BaryonField[field][i] = 0.0;
    }
  }

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  float MagneticUnits = sqrt(4.0*Pi*DensityUnits)*VelocityUnits;
  
  /* Problem parameters */
  float rho = 1.0;
  float cs = 1e-3;
  float pres = rho*cs*cs;
  float Bnaught = 0.0;
  const float q = VelocityGradient;
  const float Omega = AngularVelocity;

  /* Set up the background */

  int n = 0;  
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++, n++) {

	FLOAT x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	FLOAT y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	FLOAT z = (GridRank > 2) ? CellLeftEdge[2][k] + 0.5*CellWidth[2][k] : 0.0;
	
	BaryonField[iden][n] = rho;

	float eint, h, cs_temp, dpdrho, dpde;
	EOS(pres, BaryonField[iden][n], eint, h, cs_temp, dpdrho, dpde, EOSType, 1);
  
	float vx = 0.0;
	float vy = -x*q*Omega;
	float vz = 0.0;

	BaryonField[ivx][n] = vx;
	BaryonField[ivy][n] = vy;
	BaryonField[ivz][n] = vz;

	BaryonField[ietot][n] = eint + 0.5*(vx*vx + vy*vy + vz*vz);

	if (DualEnergyFormalism) {
	  BaryonField[ieint][n] = eint;
	}

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBz][n] = Bnaught;
	  BaryonField[ietot][n] += 0.5 * pow(Bnaught,2) / rho;
	}	

      } // end loop over grid
    }
  }

  /* ProblemType 1: Vortex wave.
     Reference: B. M. Johnson & C. F. Gammie, ApJ, 2005, 626, 978. */

  if (ShearingBoxProblemType == 1) {
    const FLOAT Lx = DomainRightEdge[0] - DomainLeftEdge[0];
    const FLOAT Ly = DomainRightEdge[1] - DomainLeftEdge[1];
    const FLOAT kx0 = -8.0*2.0*Pi/Lx;
    const FLOAT ky = 2.0*2.0*Pi/Ly;
    const float vx0 = 1e-4; // in unit of cs
    int n = 0;  
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++, n++) {
	  int igrid = i + (j + k * GridDimension[1]) * GridDimension[0];
	  FLOAT x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  FLOAT y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  
	  float dvx = vx0 * cs * cos(kx0 * x + ky * y);
	  float dvy = -kx0 / ky * dvx;
	  float drho = rho/(cs*ky) * (-2.0*kx0*kx0*q*Omega + 2.0*(q-1.0)*Omega) * vx0 * sin(kx0*x + ky*y);
	  
	  //BaryonField[iden][igrid] += drho;
	  BaryonField[ivx ][igrid] += dvx;
	  BaryonField[ivy ][igrid] += dvy;	  
	}
      }
    }
  }

  return SUCCESS;
}





