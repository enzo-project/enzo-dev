/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID HIERARCHY INFORMATION)
/
/  written by: Robert Harkness
/  date:       July, 2003
/
/  PURPOSE:
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
void my_exit(int status);
 
// function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
 
 
 
int grid::WriteStuff(FILE *fptr, char *base_name, int grid_id)
{
 
  int dim;
 
  /* make sure quantities defined at least for 3d */
 
  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }
 
  /* 1) Save general grid class data */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
 
    fprintf(fptr, "GridRank          = %"ISYM"\n", GridRank);
 
    fprintf(fptr, "GridDimension     = ");
    WriteListOfInts(fptr, GridRank, GridDimension);
 
    fprintf(fptr, "GridStartIndex    = ");
    WriteListOfInts(fptr, GridRank, GridStartIndex);
 
    fprintf(fptr, "GridEndIndex      = ");
    WriteListOfInts(fptr, GridRank, GridEndIndex);
 
    fprintf(fptr, "GridLeftEdge      = ");
    WriteListOfFloats(fptr, GridRank, GridLeftEdge);
 
    fprintf(fptr, "GridRightEdge     = ");
    WriteListOfFloats(fptr, GridRank, GridRightEdge);
 
    fprintf(fptr, "Time              = %"GOUTSYM"\n", Time);
 
    fprintf(fptr, "SubgridsAreStatic = %"ISYM"\n", SubgridsAreStatic);
 
    fprintf(fptr, "NumberOfBaryonFields = %"ISYM"\n", NumberOfBaryonFields);
 
  }
 
  if (NumberOfBaryonFields > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
 
      fprintf(fptr, "FieldType = ");
 
      WriteListOfInts(fptr, NumberOfBaryonFields, FieldType);
 
      fprintf(fptr, "BaryonFileName = %s\n", name);
 
      fprintf(fptr, "CourantSafetyNumber    = %"FSYM"\n", CourantSafetyNumber);
      fprintf(fptr, "PPMFlatteningParameter = %"ISYM"\n", PPMFlatteningParameter);
      fprintf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", PPMDiffusionParameter);
      fprintf(fptr, "PPMSteepeningParameter = %"ISYM"\n", PPMSteepeningParameter);
 
    }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  /* 3) Save particle quantities. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfParticles   = %"ISYM"\n", NumberOfParticles);
 
  if (NumberOfParticles > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(fptr, "ParticleFileName = %s\n", name); // must be same as above
 
  } // end: if (NumberOfParticles > 0)
 
  /* 4) Save Gravity info. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (SelfGravity)
      fprintf(fptr, "GravityBoundaryType = %"ISYM"\n", GravityBoundaryType);
 
  return SUCCESS;
 
}
