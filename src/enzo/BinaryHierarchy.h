/***********************************************************************
/
/  BINARY HIERARCHY STRUCTURE
/
/  written by: Matthew Turk
/  date:       October, 2009
/  modified1:
/
/  PURPOSE: The structure in which we will hold all of our arrays containing
/  grid information to output to an HDF5 file.
/
/  This will only work with 3D, but it would not really be necessary with less
/  than 3D anyway.  Note also that we only output a subset of the information
/  that goes into the standard hierarchy -- this is designed to speed up
/  analysis.
/
/  REQUIRES: macros_and_parameters.h
/
************************************************************************/

#ifndef BINARY_HIERARCHY_DEFINED__
#define BINARY_HIERARCHY_DEFINED__

struct hierarchy_arrays
{
  int grid_count;
  int current_parent;
  int *ActiveDimensions;
  FLOAT *LeftEdges;
  FLOAT *RightEdges;
  int *Level;
  int *ParentIDs;
  int *Processor;
  int *NumberOfParticles;
};

extern hierarchy_arrays HierarchyArrays;

#endif
