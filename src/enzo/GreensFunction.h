/***********************************************************************
/
/  GREENS FUNCTION CLASS
/
/  written by: Greg Bryan
/  date:       March, 1995
/
/  modified1: Jean-Claude Passy
/  date: May, 2018
/
/  PURPOSE:
/
************************************************************************/
#ifndef GREENS_FUNCTION_H
#define GREENS_FUNCTION_H

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif /* DEFINE_STORAGE */

EXTERN int GFAlreadyInitialized;      // Used for initializing
EXTERN int GFDimensions[MAX_NUMBER_OF_GREENS_FUNCTIONS][MAX_DIMENSION];
EXTERN gravity_boundary_type GFBoundary[MAX_NUMBER_OF_GREENS_FUNCTIONS];
EXTERN float *GFFunction[MAX_NUMBER_OF_GREENS_FUNCTIONS];
EXTERN int GFInUse[MAX_NUMBER_OF_GREENS_FUNCTIONS];
EXTERN int GFFunctionSize[MAX_NUMBER_OF_GREENS_FUNCTIONS];  // size in floats
EXTERN int GFIndex;                       // index for the next GF
EXTERN int GFTotalSize;                   // total number of floats used


//! Class representing a Greens function.
class greens_function
{
 private:

  float *GreensFunctionThisGrid;          // pointer for this grid.
  int    GreensFunctionIndex;

 public:

  //! Prepare the Greens functions
  int PrepareGreensFunction(int Rank, int RealDims[], int Dims[],
                            gravity_boundary_type BoundaryType,
                            int RefinementFactor);


  //! Multiplies the given (complex k-space) field with the appropriate G.F.
  int MultiplyGreensFunction(int Rank, int RealDims[], int Dims[], float *field);


  //! Compute the acceleration field over the k-space potential by multipling
  //! with D(k).  The result goes in AccelerationField. */
  int ComputeAcceleration(int Rank, int RealDims[], int Dims[], int Direction,
                          float *Potential, float *AccelerationField,
                          float CellWidth);


  //! Release a Greens function
  void Release();
};


//! Build Greens function
int MakeGreensFunction(int Rank, int RealDims[], int Dims[], float *Function,
                       gravity_boundary_type BoundaryType,
                       int RefinementFactor);

//! S2 Particle functions.  Note: 2D is NOT!!! correct.
//! It should be a generalized hypergeometric function.
float S2Function1D(float keta2);
float S2Function2D(float keta2);
float S2Function3D(float keta2);
float FloatSin(float x); // Sin function: float version.  Sometimes, C sucks.


#endif /* GREENS_FUNCTION_H */





