/***********************************************************************
/
/  RADIATION SOURCE STRUCTURE AND ROUTINES
/
/  written by: Tom Abel & Greg Bryan
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: structure for a linked list of radiation sources
/
************************************************************************/
#ifndef __RADIATIONSOURCE_H
#define __RADIATIONSOURCE_H

struct SuperSourceEntry {
  SuperSourceEntry *ParentSource;
  SuperSourceEntry *ChildSource[MAX_LEAF]; // MAX_LEAF=2 :: binary tree
  FLOAT Position[MAX_DIMENSION];
  int LeafID;
  float ClusteringRadius;
};

struct RadiationSourceEntry  {
  RadiationSourceEntry *NextSource; // Next Link
  RadiationSourceEntry *PreviousSource; // Previous Link
  SuperSourceEntry *SuperSource;  // Associated super source
  int   GridID;                   // Associated grid ID
  int   GridLevel;                // Associated grid level
  int   Type;                     // Type allows for beaming etc.      
  float Luminosity;               // Bolometric photon number luminosity
				  // in [#/s] * TimeUnits/LengthUnits^3
  float CreationTime;             // When Source is formed in code units
  float LifeTime;                 // LifeTime of source in code units
  float RampTime;                 // Time for the source to reach full luminosity
  int   EnergyBins;
  float *Energy;                  // Energy bins
  float *SED;                     // fractional Spectral energy distribution
  FLOAT *Position;                // Position of source
  float *Orientation;             // Direction for one cone of beamed rad.
  bool  AddedEmissivity;          // flag to show that we've added
                                  // emissivity for FS solver.
};

struct SuperSourceData {
  RadiationSourceEntry *Source;
  FLOAT Position[MAX_DIMENSION];
  float Luminosity;
};

#endif
