/***********************************************************************
/
/  FLUXES STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE: Structure to hold fluxes surrounding a box.  The fluxes
/    begin at a global index (i.e. an index, demarcated by zone widths
/    at this level of refinement that spans the entire grid, starting at
/    zero at DomainLeftEdge.  The left and right refer to the two faces
/    that must be recorded.  If the length of a dimension is one, then
/    no flux information is required at either end (so no space is
/    allocated for Left/RightFluxes).  The first dimension in the indecies
/    refers to the particular flux data at either end of the box in that
/    dimension (i.e. it is the same index as the second dimension of the
/    flux data themselves: Left/RightFluxes).  The second dimension
/    represents a vector that specifies the position of the first corner
/    (StartGlobalIndex) or the end corner (EndGlobalIndex).
/
/  REQUIRES: macros_and_parameters.h
/
************************************************************************/

#ifndef FLUXES_DEFINED__
#define FLUXES_DEFINED__

struct fluxes
{
  long_int LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  float *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
  float *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
};


#endif
