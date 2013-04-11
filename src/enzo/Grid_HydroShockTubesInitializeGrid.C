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

int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int grid::HydroShockTubesInitializeGrid(float x0, 
					float rhol, float rhor,
					float vxl,  float vxr,
					float vyl,  float vyr,
					float vzl,  float vzr,
					float pl,   float pr
					)
{  

  int MachNum, PSTempNum, PSDenNum;

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }

  if(ShockMethod){
    FieldType[MachNum   = NumberOfBaryonFields++] = Mach;
    if(StorePreShockFields){
      FieldType[PSTempNum = NumberOfBaryonFields++] = PreShockTemperature;
      FieldType[PSDenNum = NumberOfBaryonFields++] = PreShockDensity;
    }
  }    

  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }


  int size = 1, index, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  
  /* transform pressure to total energy */
  float etotl, etotr, v2;
  v2 = vxl * vxl + vyl * vyl + vzl * vzl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2;

  v2 = vxr * vxr + vyr * vyr + vzr * vzr;
  etotr = pr / ((Gamma-1.0)*rhor) + 0.5*v2;

  FLOAT x;
  int i;
  for (int k = 0; k < GridDimension[2]; k++) {
  for (int j = 0; j < GridDimension[1]; j++) {
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    index = GRIDINDEX_NOGHOST(i,j,k); 

    if (x <= x0) {
      BaryonField[iden ][index] = rhol;
      BaryonField[ivx  ][index] = vxl;
      BaryonField[ivy  ][index] = vyl;
      BaryonField[ivz  ][index] = vzl;
      BaryonField[ietot][index] = etotl;
      if (DualEnergyFormalism) {
	BaryonField[ieint][index] = etotl - 0.5*(vxl*vxl+vyl*vyl+vzl*vzl);
      }
    } else {
      BaryonField[iden ][index] = rhor;
      BaryonField[ivx  ][index] = vxr;
      BaryonField[ivy  ][index] = vyr;
      BaryonField[ivz  ][index] = vzr;
      BaryonField[ietot][index] = etotr;
      if (DualEnergyFormalism) {
	BaryonField[ieint][index] = etotr - 0.5*(vxr*vxr+vyr*vyr+vzr*vzr);
      }
    }

    //Shock
    if (ShockMethod) {
      BaryonField[MachNum][index] = tiny_number;
      if (StorePreShockFields) {
        BaryonField[PSTempNum][index] = tiny_number;
        BaryonField[PSDenNum][index] = tiny_number;
      }
    }

  }
  }
  }

  return SUCCESS;
}

/* Version to specify three regions */

int grid::HydroShockTubesInitializeGrid(float x0, float x1,
					float rhol, float rhor, float rhoc,
					float vxl,  float vxr, float vxc,
					float vyl,  float vyr, float vyc,
					float vzl,  float vzr, float vzc,
					float pl,   float pr, float pc
					)
{  

  int MachNum, PSTempNum, PSDenNum;

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }

  if(ShockMethod){
    FieldType[MachNum   = NumberOfBaryonFields++] = Mach;
    if(StorePreShockFields){
      FieldType[PSTempNum = NumberOfBaryonFields++] = PreShockTemperature;
      FieldType[PSDenNum = NumberOfBaryonFields++] = PreShockDensity;
    }
  }    
  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }


  int size = 1, dim, index;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];
  
  /* transform pressure to total energy */
  float etotl, etotr, etotc, v2;
  v2 = vxl * vxl + vyl * vyl + vzl * vzl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2;

  v2 = vxr * vxr + vyr * vyr + vzr * vzr;
  etotr = pr / ((Gamma-1.0)*rhor) + 0.5*v2;

  v2 = vxc * vxc + vyc * vyc + vzc * vzc;
  etotc = pc / ((Gamma-1.0)*rhoc) + 0.5*v2;

  FLOAT x;
  int i;
  for (int k = 0; k < GridDimension[2]; k++) {
  for (int j = 0; j < GridDimension[1]; j++) {
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    index = GRIDINDEX_NOGHOST(i,j,k);

    if (x <= x0) {
      BaryonField[iden ][index] = rhol;
      BaryonField[ivx  ][index] = vxl;
      BaryonField[ivy  ][index] = vyl;
      BaryonField[ivz  ][index] = vzl;
      BaryonField[ietot][index] = etotl;
      if (DualEnergyFormalism) {
	BaryonField[ieint][index] = etotl - 0.5*(vxl*vxl+vyl*vyl+vzl*vzl);
      }
    } else if (x <= x1) {
      BaryonField[iden ][index] = rhoc;
      BaryonField[ivx  ][index] = vxc;
      BaryonField[ivy  ][index] = vyc;
      BaryonField[ivz  ][index] = vzc;
      BaryonField[ietot][index] = etotc;
      if (DualEnergyFormalism) {
	BaryonField[ieint][index] = etotc - 0.5*(vxc*vxc+vyc*vyc+vzc*vzc);
      }
    }
    else {
      BaryonField[iden ][index] = rhor;
      BaryonField[ivx  ][index] = vxr;
      BaryonField[ivy  ][index] = vyr;
      BaryonField[ivz  ][index] = vzr;
      BaryonField[ietot][index] = etotr;
      if (DualEnergyFormalism) {
	BaryonField[ieint][index] = etotr - 0.5*(vxr*vxr+vyr*vyr+vzr*vzr);
      }
    }

  //Shock
    if (ShockMethod) {
      BaryonField[MachNum][index] = tiny_number;
      if (StorePreShockFields) {
	BaryonField[PSTempNum][index] = tiny_number;
	BaryonField[PSDenNum][index] = tiny_number;
      }
    }
  }
  }
  }

  return SUCCESS;
}
