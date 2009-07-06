/***********************************************************************
/
/  GRID CLASS (ACCESS THE BARYON FIELDS)
/
/  written by: Daniel R. Reynolds
/  date:       October, 2006
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
// Density field 
float* grid::AccessDensity() {
  int DensNum = -1;
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[DensNum];
}

// Total Energy field 
float* grid::AccessTotalEnergy() {
  int EnergyNum = -1;
  if ((EnergyNum = FindField(TotalEnergy,FieldType,NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[EnergyNum];
}

// Gas Energy field 
float* grid::AccessGasEnergy() {
  int GENum = -1;
  if (((GENum = FindField(InternalEnergy,FieldType,NumberOfBaryonFields))<0) 
      || (DualEnergyFormalism == FALSE)) 
    return NULL;
  return BaryonField[GENum];
}

// Velocity1 field 
float* grid::AccessVelocity1() {
  int Vel1Num = -1;
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields))<0)
    return NULL;
  return BaryonField[Vel1Num];
}

// Velocity2 field 
float* grid::AccessVelocity2() {
  int Vel2Num = -1;
  if ((Vel2Num = FindField(Velocity2, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Vel2Num];
}

// Velocity3 field 
float* grid::AccessVelocity3() {
  int Vel3Num = -1;
  if ((Vel3Num = FindField(Velocity3, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Vel3Num];
}

// Electron Density field 
float* grid::AccessElectronDensity() {
  int ENum = -1;
  if ((ENum = FindField(ElectronDensity,FieldType,NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[ENum];
}

// Hydrogen-I Density field 
float* grid::AccessHIDensity() {
  int HINum = -1;
  if ((HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HINum];
}

// Hydrogen-II Density field 
float* grid::AccessHIIDensity() {
  int HIINum = -1;
  if ((HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HIINum];
}

// Helium-I Density field 
float* grid::AccessHeIDensity() {
  int HeINum = -1;
  if ((HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeINum];
}

// Helium-II Density field 
float* grid::AccessHeIIDensity() {
  int HeIINum = -1;
  if ((HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeIINum];
}

// Helium-III Density field 
float* grid::AccessHeIIIDensity() {
  int HeIIINum = -1;
  if ((HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HeIIINum];
}

// Emissivity Field
float* grid::AccessEmissivityField()
 {
  int EtaNum = -1;
  if ((EtaNum = FindField(EmissivityField, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[EtaNum];
}

