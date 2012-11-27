/***********************************************************************
/
/  GRID CLASS (PREPARES GRID FOR INITIALIZATION)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:  Sets some basic grid quantities based on the arguments.
/
/  INPUTS:
/    GridDim   - int array of dimensions (includes ghost zones unless dim = 1)
/    LeftEdge  - float array of left edge starting positions in terms of
/                problem space (i.e. starting 'vector').  Doesn't include
/                ghost zones.
/    RightEdge - stop (right) 'vector'.
/
************************************************************************/
 
//
//  Assign basic values to a grid (allocate fields)
//
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
void grid::PrepareGrid(int Rank, int GridDim[],
		       FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles)
{
 
  /* debugging code */
/*
  if (debug) {
    printf("PrepareGrid: Rank = %"ISYM"\n", Rank);
    printf("PrepareGrid: GridDim = ");
    WriteListOfInts(stdout, Rank, GridDim);
    printf("PrepareGrid: LeftEdge = ");
    WriteListOfFloats(stdout, Rank, LeftEdge);
    printf("PrepareGrid: RightEdge = ");
    WriteListOfFloats(stdout, Rank, RightEdge);
  }
*/
  /* Set Particle quantities. */
 
  NumberOfParticles = NumParticles;
 
  /* Set global grid quantities.
     (Start/End index are zero based). */
 
  GridRank = Rank;
 
  for (int dim = 0; dim < GridRank; dim++) {
    GridDimension[dim]  = GridDim[dim];
    GridStartIndex[dim] = min(NumberOfGhostZones, GridDim[dim]-1);
    GridEndIndex[dim]   = min(ABS(GridDim[dim]-NumberOfGhostZones-1),
			      GridDim[dim]-1);
    GridLeftEdge[dim]   = LeftEdge[dim];
    GridRightEdge[dim]  = RightEdge[dim];
  }
 
  /* compute derived quantites */
 
  this->PrepareGridDerivedQuantities();

  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2 * NumberOfGhostZones);
  }
  
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }
  
  if (HydroMethod == MHD_RK) {

    if (divB == NULL) {
      divB = new float[activesize];
    }

    for (int dim = 0; dim < 3; dim++) {
      if (gradPhi[dim] == NULL) {
	gradPhi[dim] = new float[activesize];
      }
    }

    for (int dim = GridRank; dim < 3; dim++) {
      for (int n = 0; n < activesize; n++) {
	gradPhi[dim][n] = 0.0;
	divB[n] = 0.0;
      }
    }
  }


#ifdef TRANSFER
  SubgridMarker = NULL;
  if (PhotonPackages == NULL) {
    PhotonPackages = new PhotonPackageEntry;
    PhotonPackages->NextPackage = NULL;
    PhotonPackages->PreviousPackage = NULL;
  }
#endif
 
}
 
