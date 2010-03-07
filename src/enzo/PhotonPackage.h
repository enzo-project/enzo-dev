/***********************************************************************
/
/  PHOTON PACKAGE CLASS
/
/  written by: Tom Abel & Greg Bryan
/  date:       August, 2003
/  modified1:  February, 2010 by JHW
/                Converted into a poor man's class with everything 
/                public.  I need a constructor/destructor to use the
/                MemoryPool to avoid memory fragmentation.
/
/  PURPOSE: Constructs a Linked List of Photon Packages including data
/
************************************************************************/
#include "RadiationSource.h"

class PhotonPackageEntry
{
public:
  PhotonPackageEntry *NextPackage; // Next Link
  PhotonPackageEntry *PreviousPackage; // Previous Link
  SuperSourceEntry *CurrentSource;  // Currently used (super)source  
  float Photons;                // number of photons in package
  int   Type;                   // 0 = HI, 1=HeI, 2=HeII, 3=H2I_LW, 4=Xray
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

  /* CONSTRUCTOR AND DESTRUCTOR */

  PhotonPackageEntry(void);

  virtual ~PhotonPackageEntry(void) {};

  /* Overloaded new/delete to use the memory pool, if requested */

#ifdef MEMORY_POOL
  void* operator new(size_t nobjects);
  void operator delete(void* object);
#endif

};
