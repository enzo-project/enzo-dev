#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "GreensFunction.h"


/* Set the number of interpolation points per two Pi stretch in the domain. */
#define NUMBER_OF_POINTS_PER_TWOPI 20000

#define DONT_USE_LOCK
#define UNLOCKED 0
#define LOCKED 1

static int GFLock = UNLOCKED;


/* function prototypes */
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int FastFourierTransform(float *buffer, int Rank, int DimensionReal[],
             int Dimension[], int direction, int type);
int ComputeTable(float Min, float Max, float Step, float (*Function)(float),
         float **Table, float *TableStep, float *TableMin,
         float *TableMax);
extern "C" void FORTRAN_NAME(make_green)(float *green,
                                            int *nx, int *ny, int *nz,
                                            int *in, int *jn, int *kn,
                                            int *nalias,
                                         float *S2Part, int *refinement,
                                         float *S2Table, float *S2Min,
                                            float *S2Step,
                                         float *SinTable, float *SinMin,
                                            float *SinStep);

/* Definitions of the greens_function class members */
int greens_function::PrepareGreensFunction(int Rank, int RealDims[], int Dims[],
                                           gravity_boundary_type BoundaryType,
                                           int RefinementFactor)
{

#ifdef USE_LOCK
  float a;
  printf("entering GF_Prep:\n");
  while (GFLock == LOCKED) {
    printf("+");
    for (int l = 0; l < 10000; l++)
      a += pow(-1, l)*sin(l);
  }

  GFLock = LOCKED;
#endif /* USE_LOCK */

  int dim, i = -1, Match = FALSE;

  while (!Match && ++i < GreensFunctionMaxNumber) {
    if (GFBoundary[i] == BoundaryType) {
      Match = TRUE;
      for (dim = 0; dim < Rank; dim++)
    if (GFDimensions[i][dim] != Dims[dim])
      Match = FALSE;
    }
  }

  /* If found, we don't have to do any work. */

  if (Match) {
    GFInUse[i]++;
    GreensFunctionThisGrid = GFFunction[i];
    GreensFunctionIndex    = i;
    GFLock = UNLOCKED;
    return SUCCESS;
  }

  /* Error check. */

  if (GreensFunctionMaxNumber < 1) {
    fprintf(stderr, "GreensFunctionMaxNumber must be > 1!\n");
    return FAIL;
  }

  /* If not, we'll compute a new one and put it in the spot specified by
     Index (wrapping the counter if necessary). */

  if (GFIndex > GreensFunctionMaxNumber-1)
    GFIndex = min(1, GreensFunctionMaxNumber-1);  // set to 1 if possible
  int nchecks = 0;
  while (GFInUse[GFIndex] > 0 && nchecks++ < GreensFunctionMaxNumber)
    GFIndex = (GFIndex+1) % GreensFunctionMaxNumber;
  if (nchecks == GreensFunctionMaxNumber) {
    fprintf(stderr, "GF: no empty GreensFunction.\n");
    return FAIL;
  }

  /* Add one to index and set boundary to undefined so nobody uses or over-
     runs this spot. */

  int ThisGF         = GFIndex++;
  GFBoundary[ThisGF] = GravityUndefined;
  GFInUse[ThisGF]    = 1;

  /* If the spot is already being used, delete the old stuff. */

  if (GFFunction[ThisGF] != NULL) {
    delete GFFunction[ThisGF];
    GFTotalSize -= GFFunctionSize[ThisGF];   // reduce size by this amount
    if (debug) {
      printf("GreensFunction: deleting GF[%d] with dims ", ThisGF);
      WriteListOfInts(stdout, Rank, GFDimensions[ThisGF]);
    }
  }

  /* Fill in the new spot. */

  int size = RealDims[0]/2 + 1;
  for (dim = 0; dim < Rank; dim++) {
    GFDimensions[ThisGF][dim] = Dims[dim];
    if (dim > 0)
      size *= RealDims[dim];
  }
  GFFunctionSize[ThisGF] = size;
  GFTotalSize           += size;
  GFFunction[ThisGF]     = new float[size];  // allocate space for G.F.
  if (debug) {
    printf("GreensFunction: creating GF[%d] with dims ", ThisGF, size);
    WriteListOfInts(stdout, Rank, Dims);
  }
  if (GFFunction[ThisGF] == NULL) {
    fprintf(stderr, "malloc error (out of memory?)\n");
    return FAIL;
  }

  /* Create the field */

  if (BoundaryType == TopGridIsolated) {
    fprintf(stderr, "isolated gravity is not supported.\n");
    return FAIL;
    //    if (MakeGreensFunctionTopGridIsolated(Rank, RealDims, Dims,
    //			     GFFunction[ThisGF]) == FAIL) {
    //      fprintf(stderr, "Error in MakeGreensFunctionTopGridIsolated.\n");
    //      return FAIL;
    //    }
  }
  else
    if (MakeGreensFunction(Rank, RealDims, Dims, GFFunction[ThisGF],
               BoundaryType, RefinementFactor) == FAIL) {
      fprintf(stderr, "Error in MakeGreensFunction.\n");
      return FAIL;
    }

  /* Set the pointer to the new function. */

  GreensFunctionThisGrid = GFFunction[ThisGF];
  GreensFunctionIndex    = ThisGF;

  /* Set the boundary type only when done. */

  GFBoundary[ThisGF]     = BoundaryType;

  if (debug)
    printf("GreensFunction: finished creating GF[%d].\n", ThisGF);

#ifdef USE_LOCK
  GFLock = UNLOCKED;
#endif /* USE_LOCK */
  return SUCCESS;
}


int greens_function::MultiplyGreensFunction(int Rank, int RealDims[],
                                            int Dims[], float *field)
{

  int i, j, k, dim, fieldindex, greenindex;

  /* Find the appropriate Greens' function.  If it's not present, that's
     an error! */

  if (GreensFunctionThisGrid == NULL) {
    fprintf(stderr, "Green's Function not found!\n");
    return FAIL;
  }

  /* Copy Dims to FullDims and set unused indexes. */

  int FullDims[MAX_DIMENSION], FullRealDims[MAX_DIMENSION];

  for (dim = 0; dim < Rank; dim++) {
    FullDims[dim]     = Dims[dim];
    FullRealDims[dim] = RealDims[dim];
  }
  for (dim = Rank; dim < MAX_DIMENSION; dim++) {
    FullDims[dim]     = 1;
    FullRealDims[dim] = 1;
  }

  int nxmid = FullDims[0]/2 + 1;

  /* Multiply the fields. */

  for (k = 0; k < FullDims[2]; k++)
    for (j = 0; j < FullDims[1]; j++) {
      fieldindex = (k*FullRealDims[1] + j)*FullRealDims[0];
      greenindex = (k*FullRealDims[1] + j)*nxmid;
      for (i = 0; i < FullDims[0]/2+1; i++, greenindex++) {

    field[fieldindex++] *= GreensFunctionThisGrid[greenindex];
    field[fieldindex++] *= GreensFunctionThisGrid[greenindex];

      }
    } // end of loop over grids

  return SUCCESS;
}


int greens_function::ComputeAcceleration(int Rank, int RealDims[], int Dims[],
                                         int Direction, float *Potential,
                                         float *AccelerationField,
                                         float CellWidth)
{

  /* Copy Dims to FullDims and set unused indexes. */

  int dim, FullDims[MAX_DIMENSION], FullRealDims[MAX_DIMENSION];
  const float Pi = 3.14159;

  for (dim = 0; dim < Rank; dim++) {
    FullDims[dim]     = Dims[dim];
    FullRealDims[dim] = RealDims[dim];
  }
  for (dim = Rank; dim < MAX_DIMENSION; dim++) {
    FullDims[dim]     = 1;
    FullRealDims[dim] = 1;
  }

  /* Error check. */

  if (FullDims[0]+2 != FullRealDims[0] || FullDims[1] != FullRealDims[1]) {
    fprintf(stderr, "FullDims assumption failed!\n");
    return FAIL;
  }

  /* Divide by CellWidth to get the acceleration units right. */

  float factor = 2.0*Pi/float(FullDims[Direction])/CellWidth *
                 GravityResolution*GravityResolution*GravityResolution;

#ifdef USE_FORTRAN

  /* Call fortran routine to do the real work. */

  FORTRAN_NAME(mult_pot)(FullDims, FullDims+1, FullDims+2, &Direction,
                         AccelerationField, Potential, &factor);

#else /* USE_FORTRAN */

  int i, j, k, n, index = 0;
  float kvalue;

  /* (i) Direction. */

  if (Direction == 0)
    for (k = 0; k < FullDims[2]; k++)
      for (j = 0; j < FullDims[1]; j++)
    for (i = 0; i < FullDims[0]/2+1; i++, index += 2) {

      /* Compute k. */

      kvalue = float(i)*factor;

      /* Multiply potential by -ik. */

      AccelerationField[index  ] =  Potential[index+1] * kvalue;
      AccelerationField[index+1] = -Potential[index  ] * kvalue;

    }

  /* (j/k) Direction. */

  if (Direction == 1 || Direction == 2)
    for (k = 0; k < FullDims[2]; k++)
      for (j = 0; j < FullDims[1]; j++) {

    /* Compute k. */

    if (Direction == 1) n = j;
    if (Direction == 2) n = k;
    if (n > FullDims[Direction]/2)
      n -= FullDims[Direction];
    kvalue = float(n)*factor;

    /* Multiply potential by -ik. */

    for (i = 0; i < FullDims[0]/2+1; i++, index += 2) {

      AccelerationField[index  ] =  Potential[index+1] * kvalue;
      AccelerationField[index+1] = -Potential[index  ] * kvalue;

    }

      }

#endif /* USE_FORTRAN */

  return SUCCESS;
}


void greens_function::Release()
{
  GFInUse[GreensFunctionIndex]--;
};


/* Free functions and parameters definitions */
/* S2 Particle Table values. */
static float S2TableStep = 0.0;
static float *S2Table    = NULL;
static float S2TableMin  = 0.0;
static float S2TableMax  = 0.0;

/* Sin Table values. */
static float SinTableStep = 0.0;
static float *SinTable    = NULL;
static float SinTableMin  = 0.0;
static float SinTableMax  = 0.0;


int MakeGreensFunction(int Rank, int RealDims[], int Dims[], float *Function,
              gravity_boundary_type BoundaryType, int RefinementFactor)
{

#ifdef USE_LOCK
  float a;
  while (GFLock == LOCKED) {
    if (debug) printf("+");
    for (int l = 0; l < 10000; l++)
      a += pow(-1, l)*sin(l);
  }

  GFLock = LOCKED;
#endif /* USE_LOCK */

  /* declarations */

  int dim, FullDims[MAX_DIMENSION], FullRealDims[MAX_DIMENSION];

  /* Error check. */

/*  if (Rank == 2) {
    fprintf(stderr, "MakeGreen: Rank == 2 not supported.\n");
    return FAIL;
  } */

  if (BoundaryType == TopGridIsolated) {
    fprintf(stderr, "BoundaryType %d not supported.\n", BoundaryType);
    return FAIL;
  }

  if (BoundaryType == TopGridPeriodic && RefinementFactor != 1) {
    fprintf(stderr, "RefinementFactor %d incompatable with TopGridPeriodic.\n",
            RefinementFactor);
    return FAIL;
  }

  if (BoundaryType == SubGridIsolated && RefinementFactor <= 1) {
    fprintf(stderr, "RefinementFactor %d incompatable with SubgridIsolated.\n",
            RefinementFactor);
    return FAIL;
  }

  /* Copy Dims to FullDims and set unused indexes. */

  for (dim = 0; dim < Rank; dim++) {
    FullDims[dim]     = Dims[dim];
    FullRealDims[dim] = RealDims[dim];
  }
  for (dim = Rank; dim < MAX_DIMENSION; dim++) {
    FullDims[dim]     = 1;
    FullRealDims[dim] = 1;
  }

  /* Set some constants. */

  const float Pi   = 3.14159;
  int nalias = NUMBER_IN_ALIAS_SUM;

  /* Precompute the S2 particle shape look up table. */

  float MinimumStep = 2*Pi/NUMBER_OF_POINTS_PER_TWOPI;
  float MinValue = MinimumStep*0.1;
  float MaxValue = (sqrt((float(nalias)+0.5)*(float(nalias)+0.5) +
                 (float(nalias)+0.5)*(float(nalias)+0.5) +
                 (float(nalias)+0.5)*(float(nalias)+0.5) ) + 0.1)
                   *2*Pi*S2ParticleSize*0.5*float(RefinementFactor);

  float (*S2Function)(float);
  if (Rank == 1) S2Function = &S2Function1D;
  if (Rank == 2) S2Function = &S2Function2D;
  if (Rank == 3) S2Function = &S2Function3D;

  if (ComputeTable(MinValue, MaxValue, MinimumStep, S2Function,
           &S2Table, &S2TableStep, &S2TableMin, &S2TableMax)
      == FAIL) {
    fprintf(stderr, "Error in ComputeTable (S2Table).\n");
    return FAIL;
  }

  /* Precompute the sin function look up table. */

  MinimumStep = 2*Pi/NUMBER_OF_POINTS_PER_TWOPI;
  MinValue = -(2.1 + float(nalias))*Pi;
  MaxValue =  (2.1 + float(nalias))*Pi;

  if (ComputeTable(MinValue, MaxValue, MinimumStep, &FloatSin,
           &SinTable, &SinTableStep, &SinTableMin, &SinTableMax)
      == FAIL) {
    fprintf(stderr, "Error in ComputeTable (SinTable).\n");
    return FAIL;
  }

  /* Call fortran function to compute Green's function. */

  FORTRAN_NAME(make_green)(Function, FullDims, FullDims+1, FullDims+2,
                           FullRealDims, FullRealDims+1, FullRealDims+2,
                           &nalias, &S2ParticleSize, &RefinementFactor,
                           S2Table, &S2TableMin, &S2TableStep,
                           SinTable, &SinTableMin, &SinTableStep);

/* #define UNUSED */
#ifdef UNUSED

  /* Check to see what this looks like in the real domain. */

  int i, j, k, size = 1, nxmid = FullDims[0]/2 + 1;
  for (dim = 0; dim < Rank; dim++)
    size *= FullRealDims[dim];

  float *test = new float[size];

  for (k = 0; k < FullDims[2]; k++)
    for (j = 0; j < FullDims[1]; j++)
      for (i = 0; i < nxmid; i++) {
    *(test + k*FullRealDims[0]*FullRealDims[1]
           + j*FullRealDims[0] + i*2) =
    *(Function + k*nxmid*FullRealDims[1] + j*nxmid + i);
    *(test + k*FullRealDims[0]*FullRealDims[1]
           + j*FullRealDims[0] + i*2+1) = 0;
      }

  FastFourierTransform(test, Rank, FullRealDims, FullDims, FFT_INVERSE, REAL_TO_COMPLEX);

  if (BoundaryType == SubGridIsolated) {
    FILE *fptr = fopen("Green.out","w");
    for (k = 0; k < FullDims[2]; k++)
      for (j = 0; j < FullDims[1]; j++)
    for (i = 0; i < FullDims[0]; i++)
      fprintf(fptr, "%e\n", *(test + k*FullRealDims[0]*FullRealDims[1]
                           + j*FullRealDims[0] + i));
    fclose(fptr);
  }

  delete test;

#endif /* UNUSED */

#ifdef USE_LOCK
  GFLock = UNLOCKED;
#endif /* USE_LOCK */

  return SUCCESS;
}


float S2Function1D(float keta2)
{
  return 2.0 / (keta2*keta2) * (1 - cos(keta2));
}


float S2Function2D(float keta2)
{
  return 12.0 / pow(keta2, 4) * (2 - 2*cos(keta2) - keta2*sin(keta2));
}


float S2Function3D(float keta2)
{
  return 12.0 / pow(keta2, 4) * (2 - 2*cos(keta2) - keta2*sin(keta2));
}


float FloatSin(float x)
{
  return float(sin(double(x)));
}
