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

  int DensNum = -1;
  int Vel1Num = -1;
  int Vel2Num = -1;
  int Vel3Num = -1;
  int TENum = -1;
  int IENum = -1;
  int MachNum = -1;
  int PSTempNum = -1;
  int PSDenNum = -1;

  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[Vel1Num = NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2) {
    FieldType[Vel2Num = NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2 || HydroMethod > 2) {
      FieldType[Vel3Num = NumberOfBaryonFields++] = Velocity3;
    }
  }
  FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[IENum = NumberOfBaryonFields++] = InternalEnergy;
  }

  if(ShockMethod){
    FieldType[MachNum = NumberOfBaryonFields++] = Mach;
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
  this->AllocateGrids();

  
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
      BaryonField[DensNum][index] = rhol;
      BaryonField[Vel1Num][index] = vxl;
      if (Vel2Num > -1) {
        BaryonField[Vel2Num][index] = vyl;
      }
      if (Vel3Num > -1) {
        BaryonField[Vel3Num][index] = vzl;
      }
      BaryonField[TENum][index] = etotl;
      if (DualEnergyFormalism) {
	BaryonField[IENum][index] = etotl - 0.5*(vxl*vxl+vyl*vyl+vzl*vzl);
      }
    } else {
      BaryonField[DensNum][index] = rhor;
      BaryonField[Vel1Num][index] = vxr;
      if (Vel2Num > -1) {
        BaryonField[Vel2Num][index] = vyr;
      }
      if (Vel3Num > -1) {
        BaryonField[Vel3Num][index] = vzr;
      }
      BaryonField[TENum][index] = etotr;
      if (DualEnergyFormalism) {
	BaryonField[IENum][index] = etotr - 0.5*(vxr*vxr+vyr*vyr+vzr*vzr);
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

  int DensNum = -1;
  int Vel1Num = -1;
  int Vel2Num = -1;
  int Vel3Num = -1;
  int TENum = -1;
  int IENum = -1;
  int MachNum = -1;
  int PSTempNum = -1;
  int PSDenNum = -1;

  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[Vel1Num = NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || HydroMethod > 2) {
    FieldType[Vel2Num = NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2 || HydroMethod > 2) {
      FieldType[Vel3Num = NumberOfBaryonFields++] = Velocity3;
    }
  }
  FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[IENum = NumberOfBaryonFields++] = InternalEnergy;
  }

  if (ShockMethod) {
    FieldType[MachNum = NumberOfBaryonFields++] = Mach;
    if (StorePreShockFields) {
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
      BaryonField[DensNum][index] = rhol;
      BaryonField[Vel1Num][index] = vxl;
      if (Vel2Num > -1) {
        BaryonField[Vel2Num][index] = vyl;
      }
      if (Vel3Num > -1) {
        BaryonField[Vel3Num][index] = vzl;
      }
      BaryonField[TENum][index] = etotl;
      if (DualEnergyFormalism) {
	BaryonField[IENum][index] = etotl - 0.5*(vxl*vxl+vyl*vyl+vzl*vzl);
      }
    } else if (x <= x1) {
      BaryonField[DensNum][index] = rhoc;
      BaryonField[Vel1Num][index] = vxc;
      if (Vel2Num > -1) {
        BaryonField[Vel2Num][index] = vyc;
      }
      if (Vel3Num > -1) {
        BaryonField[Vel3Num][index] = vzc;
      }
      BaryonField[TENum][index] = etotc;
      if (DualEnergyFormalism) {
	BaryonField[IENum][index] = etotc - 0.5*(vxc*vxc+vyc*vyc+vzc*vzc);
      }
    } else {
      BaryonField[DensNum][index] = rhor;
      BaryonField[Vel1Num][index] = vxr;
      if (Vel2Num > -1) {
        BaryonField[Vel2Num][index] = vyr;
      }
      if (Vel3Num > -1) {
        BaryonField[Vel3Num][index] = vzr;
      }
      BaryonField[TENum][index] = etotr;
      if (DualEnergyFormalism) {
	BaryonField[IENum][index] = etotr - 0.5*(vxr*vxr+vyr*vyr+vzr*vzr);
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
