/***********************************************************************
/
/  GRID CLASS (INITIALIZE MHD 3D TEST)
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
int grid::MHD3DTestInitializeGrid(int MHD3DProblemType,
				  float rhol, float rhou,
				  float vxl,  float vxu,
				  float vyl,  float vyu,
				  float pl,   float pu,
				  float Bxl,  float Bxu,
				  float Byl,  float Byu)
{  

  /* create fields */
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

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1, activesize = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);
  
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
    }
  }

  /* transform pressure to total energy */
  float etotl, etotu, v2, B2;
  v2 = vxl * vxl + vyl * vyl;
  B2 = Bxl * Bxl + Byl * Byl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2 + 0.5*B2/rhol;

  v2 = vxu * vxu + vyu * vyu;
  B2 = Bxu * Bxu + Byu * Byu;
  etotu = pu / ((Gamma-1.0)*rhou) + 0.5*v2 + 0.5*B2/rhou;

  
  int igrid;
  FLOAT x, y, z;
  if (MHD3DProblemType == 0) { // Planar shock
  float pres, eintl, eintu, h, cs, dpdrho, dpde;
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
      
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (x <= 0.1) {
	  if (k==2) {
	    printf("rhol=%"GSYM", vxl=%"GSYM", pl=%"GSYM"\n", rhol, vxl, pl);
	  }

	  BaryonField[iden ][igrid] = rhol;
	  BaryonField[ivx  ][igrid] = vxl;
	  BaryonField[ivy  ][igrid] = vyl;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotl;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pl / ((Gamma-1.0)*rhol);
	  }
	  BaryonField[iBx  ][igrid] = Bxl;
	  BaryonField[iBy  ][igrid] = Byl;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	} else {
	  if (k==2) {
	    printf("rhor=%"GSYM", vxr=%"GSYM", pr=%"GSYM"\n", rhou, vxu, pu);
	  }
	  BaryonField[iden ][igrid] = rhou;
	  BaryonField[ivx  ][igrid] = vxu;
	  BaryonField[ivy  ][igrid] = vyu;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotu;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pu / ((Gamma-1.0)*rhou);
	  }
	  BaryonField[iBx  ][igrid] = Bxu;
	  BaryonField[iBy  ][igrid] = Byu;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
      }
    }  
  }
  }

  if (MHD3DProblemType == 1) { // Uniform Density with a Shear
    float pres, eintl, eintu, h, cs, dpdrho, dpde;
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  /* Compute position */
	  igrid = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	  
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][j] + 0.5*CellWidth[2][j];
	  
	  float rho;
	  rho = 1.;
	  pres = rho/Gamma; // sound speed = 1
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	  // impose mode perturbation
	  vxl = 0.5*cos(2.0*M_PI*y);
	  vyl = 0.;
	  etotl = eintl + 0.5*(vxl*vxl + vyl*vyl) + 0.5*(Bxl*Bxl+Byl*Byl)/rho;
	  BaryonField[iden ][igrid] = rho;
	  BaryonField[ivx  ][igrid] = vxl;
	  BaryonField[ivy  ][igrid] = vyl;
	  BaryonField[ivz  ][igrid] = 0.0;
	  
	  BaryonField[ietot][igrid] = etotl;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pl / ((Gamma-1.0)*rho);
	  }
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx  ][igrid] = Bxl;
	    BaryonField[iBy  ][igrid] = Byl;
	    BaryonField[iBz  ][igrid] = 0.0;
	    BaryonField[iPhi ][igrid] = 0.0;
	  }
	  
	}
      }  
    }
  }
  
  
  
  
  return SUCCESS;
}
