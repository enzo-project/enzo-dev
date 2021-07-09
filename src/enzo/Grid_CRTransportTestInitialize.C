#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "phys_constants.h"
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

int grid::CRTransportTestInitializeGrid(int test_type, float center, 
					float rho, float vx,
					float vy,  float vz,
					float pg,  float ecr,
					float bx,  float by,
					float bz)
{
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[NumberOfBaryonFields++] = InternalEnergy;

  if (HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  FieldType[NumberOfBaryonFields++] = CRDensity;

  int iCRD = FindField(CRDensity, FieldType, NumberOfBaryonFields);

  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1, activesize = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];
  
  /* transform pressure to total energy */
  float etot, Ecr, v2, Bx, By,Bz, B2, gauss, r2, r, phi, Rho;
  v2 = vx * vx + vy * vy + vz * vz;
  etot = pg; // / ((Gamma-1.0)*rho) + 0.5*v2;

  FLOAT x, y, z;
  float D = 0.05;
  int i, j, k, igrid;
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = 0;
	z = 0;
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	 //	 r2 = (x - center) * (x - center) * (y - center)*(y - center) + (z - center) * (z - center); 
	 r2 = x*x + y*y + z*z;
	 phi = fabs(atan2(y, x));
	 r = sqrt(r2);
	 
	 // by default, set up uniform values
	 Ecr = ecr; 
	 Bx = bx;
	 By = by;
	 Bz = bz; 
	 Rho = rho;

	 if (test_type == 0){
	   Ecr = 1 + ecr*exp(-r2/(2*D));
	   Bx = bx; 
	   By = by; 
	   Bz = bz; 
	 }

	 // anisotropic diffusion
	 else if (test_type == 1){
	   if (r > 0.01){
	     Bx = -y*bx/r;
	     By = x*by/r;
	     Bz = bz; 
	   }
	   else{
	     Bx = sign(y)*bx;
	     By = by;
	     Bz = bz;
	   }
	   if (bx * by == 0){
	     Bx = bx;
	     By = by;
	   }
	   
	   if ((r < 0.7) && (r > 0.5) && (phi < pi/12.0))
	     Ecr = ecr * 1.2;
	   else
	     Ecr = ecr;
	 }
	 
	 // simple cr streaming in 1d
	 else if (test_type == 2) {
	   Ecr = 2 - fabs(x);
	   Ecr = max(Ecr, 0);
	 }
	 
	 // cr streaming bottleneck 1d
	 // described in appendix of Jiang+Oh 2018
	 else if (test_type == 3) { 
	   Rho = 0.1 + (1.0 - 0.1) * (1.0 + tanh((x - 200.)/25.)) * (1.0 + tanh((200.-x) / 25.));
	   Ecr = 1e-6;
	   if (i < GridStartIndex[0])
	     Ecr = 3.0;
	 }

	 BaryonField[iden ][igrid] = Rho;
	 BaryonField[ivx  ][igrid] = vx;
	 BaryonField[ivy  ][igrid] = vy;
	 BaryonField[ivz  ][igrid] = vz;
	 BaryonField[ietot][igrid] = etot;
	 if (CRModel)
	   BaryonField[iCRD][igrid] = Ecr;

	 if (DualEnergyFormalism) {
	   BaryonField[ieint][igrid] = pg; // / ((Gamma-1.0)*rho);
	 }
	 if (HydroMethod == MHD_RK) {
	   B2 = Bx*Bx + By*By + Bz*Bz; 
	   BaryonField[iBx][igrid] = Bx;
	   BaryonField[iBy][igrid] = By;
	   BaryonField[iBz][igrid] = Bz;
	   BaryonField[iPhi][igrid] = 0.0;
	   BaryonField[iEtot][igrid] += 0.5* B2/Rho; 
	 }
       } 
     }
  }

  return SUCCESS;
}
