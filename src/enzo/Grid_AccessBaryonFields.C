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
#include "ErrorExceptions.h"
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

// Molecular Hydrogen Density field 
float* grid::AccessHMDensity() {
  int HMNum = -1;
  if ((HMNum = FindField(HMDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HMNum];
}

// H2I Density field 
float* grid::AccessH2IDensity() {
  int H2INum = -1;
  if ((H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[H2INum];
}

// H2II Density field 
float* grid::AccessH2IIDensity() {
  int H2IINum = -1;
  if ((H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[H2IINum];
}

// DI Density field 
float* grid::AccessDIDensity() {
  int DINum = -1;
  if ((DINum = FindField(DIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[DINum];
}

// DII Density field 
float* grid::AccessDIIDensity() {
  int DIINum = -1;
  if ((DIINum = FindField(DIIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[DIINum];
}

// HDI Density field 
float* grid::AccessHDIDensity() {
  int HDINum = -1;
  if ((HDINum = FindField(HDIDensity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[HDINum];
}

// SN Colour field 
float* grid::AccessSNColour() {
  int SNColourNum = -1;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[SNColourNum];
}

// Metallicity field 
float* grid::AccessMetallicity() {
  int MetallicityNum = -1;
  if ((MetallicityNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[MetallicityNum];
}

// ExtraType0 field 
float* grid::AccessExtraType0() {
  int ExtraNum = -1;
  if ((ExtraNum = FindField(ExtraType0, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[ExtraNum];
}

// ExtraType1 field 
float* grid::AccessExtraType1() {
  int ExtraNum = -1;
  if ((ExtraNum = FindField(ExtraType1, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[ExtraNum];
}

// kphHI field 
float* grid::AccessKPhHI() {
  int kphHINum = -1;
  if ((kphHINum = FindField(kphHI, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHINum];
}

// PhotoGamma field 
float* grid::AccessPhotoGamma() {
  int PhotoGammaNum = -1;
  if ((PhotoGammaNum = FindField(PhotoGamma, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[PhotoGammaNum];
}

// kphHeI field 
float* grid::AccessKPhHeI() {
  int kphHeINum = -1;
  if ((kphHeINum = FindField(kphHeI, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHeINum];
}

// Gamma HeI field 
float* grid::AccessGammaHeI() {
  int GammaHeINum = -1;
  if ((GammaHeINum = FindField(gammaHeI, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[GammaHeINum];
}

// kphHeII field 
float* grid::AccessKPhHeII() {
  int kphHeIINum = -1;
  if ((kphHeIINum = FindField(kphHeII, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kphHeIINum];
}

// Gamma HeII field 
float* grid::AccessGammaHeII() {
  int GammaHeIINum = -1;
  if ((GammaHeIINum = FindField(gammaHeII, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[GammaHeIINum];
}

// kdissH2I field 
float* grid::AccessKDissH2I() {
  int kdissH2Inum = -1;
  if ((kdissH2Inum = FindField(kdissH2I, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[kdissH2Inum];
}

// Gravitational Potential field 
float* grid::AccessGravPotential() {
  int GravPotNum = -1;
  if ((GravPotNum = FindField(GravPotential, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[GravPotNum];
}

// Acceleration0 field 
float* grid::AccessAcceleration0() {
  int Acceleration0Num = -1;
  if ((Acceleration0Num = FindField(Acceleration0, FieldType, 
				    NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Acceleration0Num];
}

// Acceleration1 field 
float* grid::AccessAcceleration1() {
  int Acceleration1Num = -1;
  if ((Acceleration1Num = FindField(Acceleration1, FieldType, 
				    NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Acceleration1Num];
}

// Acceleration2 Density field 
float* grid::AccessAcceleration2() {
  int Acceleration2Num = -1;
  if ((Acceleration2Num = FindField(Acceleration2, FieldType, 
				    NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[Acceleration2Num];
}

// Radiation Pressure0 field 
float* grid::AccessRadPressure0() {
  int RadPressure0Num = -1;
  if ((RadPressure0Num = FindField(RadPressure0, FieldType, 
				   NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadPressure0Num];
}

// Radiation Pressure1 field 
float* grid::AccessRadPressure1() {
  int RadPressure1Num = -1;
  if ((RadPressure1Num = FindField(RadPressure1, FieldType, 
				   NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadPressure1Num];
}

// Radiation Pressure2 field 
float* grid::AccessRadPressure2() {
  int RadPressure2Num = -1;
  if ((RadPressure2Num = FindField(RadPressure2, FieldType, 
				   NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadPressure2Num];
}

// Radiation Energy (Grey, or 0th group)
float* grid::AccessRadiationFrequency0() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq0, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (1st group)
float* grid::AccessRadiationFrequency1() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq1, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (2nd group)
float* grid::AccessRadiationFrequency2() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq2, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (3rd group)
float* grid::AccessRadiationFrequency3() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq3, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (4th group)
float* grid::AccessRadiationFrequency4() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq4, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (5th group)
float* grid::AccessRadiationFrequency5() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq5, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (6th group)
float* grid::AccessRadiationFrequency6() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq6, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (7th group)
float* grid::AccessRadiationFrequency7() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq7, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (8th group)
float* grid::AccessRadiationFrequency8() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq8, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Radiation Energy (9th group)
float* grid::AccessRadiationFrequency9() {
  int RadNum = -1;
  if ((RadNum = FindField(RadiationFreq9, FieldType, 
			  NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[RadNum];
}

// Emissivity0 field 
float* grid::AccessEmissivity0() {
  int EtaNum = -1;
  if ((EtaNum = FindField(Emissivity0, FieldType, NumberOfBaryonFields))<0) 
    return NULL;
  return BaryonField[EtaNum];
}

