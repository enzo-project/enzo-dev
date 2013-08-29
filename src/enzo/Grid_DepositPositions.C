/***********************************************************************
/
/  GRID CLASS (DEPOSIT POSITIONS ONTO THE GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
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
 
extern "C" void PFORTRAN_NAME(cic_deposit)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        float *densfield, float *field, FLOAT *leftedge,
			int *dim1, int *dim2, int *dim3, float *cellsize,
					   float *cloudsize);
 
extern "C" void PFORTRAN_NAME(smooth_deposit)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        float *densfield, float *field, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, float *cellsize,
			       float *rsmooth);
 
 
int grid::DepositPositions(FLOAT *Position[], float *Mass, int Number,
			   int DepositField)
{
  if (Number == 0) return SUCCESS;

  /* DepositField specifies where the particles should go.  Set LeftEdge,
     Dimension, CellSize, DepositFieldPointer according to it. */
 
  float *DepositFieldPointer, CellSize, CloudSize;
  FLOAT LeftEdge[MAX_DIMENSION];
  int   dim, Dimension[MAX_DIMENSION];
 
  /* 1) GravitatingMassField. */
 
  if (DepositField == GRAVITATING_MASS_FIELD) {
    DepositFieldPointer = GravitatingMassField;
    CellSize            = float(GravitatingMassFieldCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = GravitatingMassFieldLeftEdge[dim];
      Dimension[dim] = GravitatingMassFieldDimension[dim];
    }
  }
 
  /* 2) GravitatingMassFieldParticles. */
 
  else if (DepositField == GRAVITATING_MASS_FIELD_PARTICLES) {
    DepositFieldPointer = GravitatingMassFieldParticles;
    CellSize            = float(GravitatingMassFieldParticlesCellSize);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = GravitatingMassFieldParticlesLeftEdge[dim];
      Dimension[dim] = GravitatingMassFieldParticlesDimension[dim];
    }
  }
 
  /* 3) MassFlaggingField */
 
  else if (DepositField == MASS_FLAGGING_FIELD) {
    DepositFieldPointer = MassFlaggingField;
    CellSize            = float(CellWidth[0][0]);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = CellLeftEdge[dim][0];
      Dimension[dim] = GridDimension[dim];
    }
  }

  /* 4) ParticleMassFlaggingField */
 
  else if (DepositField == PARTICLE_MASS_FLAGGING_FIELD) {
    DepositFieldPointer = ParticleMassFlaggingField;
    CellSize            = float(CellWidth[0][0]);
    for (dim = 0; dim < GridRank; dim++) {
      LeftEdge[dim]  = CellLeftEdge[dim][0];
      Dimension[dim] = GridDimension[dim];
    }
  }
 
  /* 5) error */
 
  else {
    ENZO_VFAIL("DepositField = %"ISYM" not recognized.\n", DepositField)
  }
 
  /* Error check. */
 
  if (DepositFieldPointer == NULL) {
    ENZO_VFAIL("DepositFieldPointer (%"ISYM") is NULL, Number = %"ISYM".\n",
	    DepositField, Number)
  }
 
  if (GridRank != 3) {
    ENZO_FAIL("New gravity module currently supports only 3d.\n");
  }
 
  if (DepositPositionsParticleSmoothRadius < CellSize)

  {
    /* Deposit to field using CIC. */
 
//  fprintf(stderr, "------DP Call Fortran cic_deposit with CellSize = %"GSYM"\n", CellSize);
    float CloudSize = CellSize;  // we assume deposit is only on self
 
    PFORTRAN_NAME(cic_deposit)(Position[0], Position[1], Position[2], &GridRank,
			      &Number, Mass, DepositFieldPointer, LeftEdge,
			      Dimension, Dimension+1, Dimension+2,
			       &CellSize, &CloudSize);
  }
  else
  {
    /* Deposit to field using large-spherical CIC, with radius of
       DepositPositionsParticleSmoothRadius */
 
//  fprintf(stderr, "------DP Call Fortran smooth_deposit with DPPSmoothRadius = %"GSYM"\n", DepositPositionsParticleSmoothRadius);
 
    PFORTRAN_NAME(smooth_deposit)(
			  Position[0], Position[1], Position[2], &GridRank,
			  &Number, Mass, DepositFieldPointer, LeftEdge,
			  Dimension, Dimension+1, Dimension+2,
			  &CellSize, &DepositPositionsParticleSmoothRadius);
  }
 
  return SUCCESS;
}
