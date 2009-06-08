/***********************************************************************
/
/  PHOTON PACKAGE STRUCTURE AND ROUTINES
/
/  written by: Tom Abel & Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: Constructs a Linked List of Photon Packages including data
/
************************************************************************/
#include "RadiationSource.h"

struct PhotonPackageEntry  {
  PhotonPackageEntry *NextPackage; // Next Link
  PhotonPackageEntry *PreviousPackage; // Previous Link
  SuperSourceEntry *CurrentSource;  // Currently used (super)source  
  float Photons;                // number of photons in package
  int   Type;                   // 0 = HI, 1=HeI, 2=HeII, 3=Xray, 4=H2I_LW
  float Energy;                 // Mean energy of photons in this package [eV]
  double CrossSection;          // Cross-section of absorber [cm^2]
  FLOAT EmissionTimeInterval;   // Time over which package was emitted
  FLOAT EmissionTime;           // Time when package was emitted
  FLOAT CurrentTime;            // Current Time
  FLOAT Radius;                 // Distance travelled
  float ColumnDensity;           // Column Density (for shielding functions)
  long  ipix;                   // pixel in HEALPIX terminology
  int   level;                  // level in HEALPIX terminology
  FLOAT SourcePosition[3];      // Position where package was emitted
  float SourcePositionDiff;     // Radius at which it was radiated (0 = pt src)
};
