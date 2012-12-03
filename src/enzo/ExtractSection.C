/***********************************************************************
/
/  EXTRACTS A SECTION OF THE GRID AND SAVES IT
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
************************************************************************/
 
// This function extracts a section of the solution at the specified level
//   and with the specified size.  The solution is taken from grids at that
//   level where possible and interpolated from above otherwise.
 
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
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
 
/* function prototypes */
 
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void my_exit(int status);
 
/* Information for creating graphical grid structure. */
 
int NumberOfLinks[3] = {2, 5, 16};
int LinkSide[16][3] = {
  {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 0},
  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1},
  {1, 0, 1}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1},
  {1, 0, 0}
};
char ExtractOutputName[] = "amr_extract";
 
 
void ExtractSection(HierarchyEntry &TopGrid, TopGridData &MetaData,
		    LevelHierarchyEntry *LevelArray[],
		    ExternalBoundary *Exterior,
		    int ExtractStart[], int ExtractEnd[],
		    FLOAT ExtractStartCoordinate[],
		    FLOAT ExtractEndCoordinate[], int ExtractLevel)
{
 
  int i, dim, level, TotalRefineBy = 1;
  int ExtractDims[MAX_DIMENSION], SetTopGridBoundary = TRUE;
  FLOAT LeftPosition[MAX_DIMENSION], RightPosition[MAX_DIMENSION],
        TempCellWidth;
  LevelHierarchyEntry *Temp;
 
  /* If undefined, set parameters. */
 
  if (ExtractLevel == INT_UNDEFINED)
    ExtractLevel = 0;
 
  for (i = 0; i < ExtractLevel; i++)
    TotalRefineBy *= RefineBy;
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
 
    /* If the start/end coordinate have been set, use them to set the
       indexes. */
 
    TempCellWidth = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
      FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);
 
    if (ExtractStartCoordinate[dim] != FLOAT_UNDEFINED)
      ExtractStart[dim] = nint((ExtractStartCoordinate[dim] -
				DomainLeftEdge[dim] ) / TempCellWidth );
 
    if (ExtractEndCoordinate[dim] != FLOAT_UNDEFINED)
      ExtractEnd[dim] = nint((ExtractEndCoordinate[dim] -
			      DomainLeftEdge[dim] ) / TempCellWidth ) - 1;
 
    /* If start/end indexes haven't been set, then set some default
       values. */
 
    if (ExtractStart[dim] == INT_UNDEFINED)
      ExtractStart[dim] = 0;
    if (ExtractEnd[dim] == INT_UNDEFINED)
      ExtractEnd[dim] = MetaData.TopGridDims[dim]*TotalRefineBy - 1;
 
    /* Unless the entire region is being extracted, we can't set the
       external boundary. */
 
    if (ExtractStart[dim] != 0 ||
	ExtractEnd[dim] != MetaData.TopGridDims[dim]*TotalRefineBy - 1)
      SetTopGridBoundary = FALSE;
 
  }
 
  /* For each grid on this level collect all the particles below it.
     Notice that this must be done even for static hierarchy's.  */
 
  for (level = MAX_DEPTH_OF_HIERARCHY-1; level > 0; level--) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridHierarchyEntry->ParentGrid->GridData->MoveAllParticles
	                                      (1, &Temp->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->MoveAllParticles.\n");
	my_exit(EXIT_FAILURE);
      }
      if (ProblemType != 60) //AK no unwanted output, please...
	printf("Called MoveAllParticles \n");
      Temp = Temp->NextGridThisLevel;
    }
  } // end: loop over levels
 
#define NO_CREATE_DENSITY_SQUARED
#ifdef CREATE_DENSITY_SQUARED
 
  /* First, modify each grid to create a new field -- rho^2. */
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->CreateDensitySquaredField();
      Temp = Temp->NextGridThisLevel;
    }
  }
 
  /* Next, project data from high levels to low levels. */
 
  for (level = MAX_DEPTH_OF_HIERARCHY-1; level > 0; level--) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ProjectSolutionToParentGrid(
	   *Temp->GridHierarchyEntry->ParentGrid->GridData);
      Temp = Temp->NextGridThisLevel;
    }
  }
 
#endif
 
  /* --------------------------------------------------------------- */
  /* Create and clear a set of grids. */
 
  grid* Grids[MAX_DEPTH_OF_HIERARCHY];
  for (i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++)
    Grids[i] = NULL;
 
  /* Loop up hierarchy, starting from ExtractLevel, creating grids until we
     get to the top. */
 
  for (level = ExtractLevel; level >= 0; level--) {
 
    /* Convert Start/End indexes into numbers and add buffer zones. */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
 
      /* Make sure the indexes are on a parental boundary. */
 
      if (level > 0) {
	ExtractStart[dim] =  int(ExtractStart[dim]/RefineBy)   *RefineBy;
	ExtractEnd[dim]   = (int(ExtractEnd[dim]  /RefineBy)+1)*RefineBy-1;
      }
 
      /* Compute corresponding positions and compute dims. */
 
      LeftPosition[dim] = FLOAT(ExtractStart[dim])/
	                  FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy)
			    *(DomainRightEdge[dim] - DomainLeftEdge[dim]) +
			      DomainLeftEdge[dim];
      RightPosition[dim] = FLOAT(ExtractEnd[dim]+1)/
                          FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy)
			    *(DomainRightEdge[dim] - DomainLeftEdge[dim]) +
			      DomainLeftEdge[dim];
      ExtractDims[dim] = ExtractEnd[dim] - ExtractStart[dim] + 1
	   + 2*NumberOfGhostZones;
    }
 
    /* Create a new grid and fill it full of goodies. */
 
    Grids[level] = new grid;
    Grids[level]->InheritProperties(LevelArray[0]->GridData);
    Grids[level]->PrepareGrid(MetaData.TopGridRank, ExtractDims,
			      LeftPosition, RightPosition, 0);
    Grids[level]->AllocateGrids();
    Grids[level]->SetTime(MetaData.Time);
 
    /* Next level up (change start/stop indicies). */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      ExtractStart[dim] = ExtractStart[dim]/RefineBy;
      ExtractEnd[dim]   = (ExtractEnd[dim]+1)/RefineBy - 1;
    }
 
    TotalRefineBy /= RefineBy;
 
  } // end loop over levels
 
  /* --------------------------------------------------------------- */
  /* Reverse direction, interpolating & copying grid values from the hierarchy
     to our set of grids. */
 
  for (level = 0; level <= ExtractLevel; level++) {
 
    /* If level > 0: Interpolate from higher grid (if present). */
 
    if (level > 0)
      if (Grids[level]->InterpolateFieldValues(Grids[level-1]) == FAIL) {
	fprintf(stderr, "Error in grid->InterpolateFieldValues.\n");
	my_exit(EXIT_FAILURE);
      }
 
    /* Copy zones from other grids on this level. */
 
    if (CopyOverlappingZones(Grids[level], &MetaData, LevelArray, level)
	== FAIL) {
      fprintf(stderr, "Error in CopyOverlappingZones.\n");
      my_exit(EXIT_FAILURE);
    }
 
    /* Set boundary values for the top grid (only for the top level and then
       only if it is the entire level). */
 
//#ifdef HENRY_FIX
    if (level == 0 && SetTopGridBoundary == TRUE)
      if (Grids[0]->SetExternalBoundaryValues(Exterior) == FAIL) {
	fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
	my_exit(EXIT_FAILURE);
      }
//#endif /* HENRY_FIX */
 
  } // end loop over levels
 
  /* Move particles from the top grid of the input hierarchy into the bottom
     grid of our fake hierarchy. */
 
// Count maximum number of particles at this level
 
  int RHC = 0;
  int RHN = 0;
 
  Temp = LevelArray[0];
 
  while (Temp != NULL) {
    RHN = Temp->GridData->ReturnNumberOfParticles();
    RHC = RHC + RHN;
    Temp = Temp->NextGridThisLevel;
  }
  printf("Allocate N %"ISYM"\n",RHC);
 
// Allocate arrays for entire particle subgrid data
 
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION];
  float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  float *Mass;
  PINT  *Number;
  int   *Type;
 
  Mass = new float[RHC];
  Number = new PINT[RHC];
  Type = new int[RHC];
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
  {
    Position[dim] = new FLOAT[RHC];
    Velocity[dim] = new float[RHC];
  }
  for (i = 0; i < NumberOfParticleAttributes; i++)
  {
    Attribute[i] = new float[RHC];
  }
 
//
 
  int RHP = 0; // global counter
 
//
 
  Temp = LevelArray[0];
 
  while (Temp != NULL) {
    if (Grids[ExtractLevel]->MoveSubgridParticles(Temp->GridData,
                                                  &RHP,
                                                  Number,
                                                  Type,
                                                  Mass,
                                                  Position,
                                                  Velocity,
                                                  Attribute) == FAIL)
    {
      fprintf(stderr, "Error in grid->MoveSubgridParticles.\n");
      my_exit(EXIT_FAILURE);
    }
    printf("Called MoveSubgridParticles %"ISYM"\n",RHP);
    Temp = Temp->NextGridThisLevel;
  }
 
//  RHP = RHP - 1;
  Grids[ExtractLevel]->SetParticlePointers(Mass, Number, Type, Position, Velocity, Attribute);
  Grids[ExtractLevel]->SetNumberOfParticles(RHP);
 
/*
  delete Mass;
  delete Number;
  delete Type;
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
  {
    delete Position[dim];
    delete Velocity[dim];
  }
  for (i = 0; i < NumberOfParticleAttributes; i++)
  {
    delete Attribute[i];
  }
*/
 
  /* Write out bottom-most grid. */
 
  FILE *output = stdout;
 
  if (debug)
    output = stderr;
 
  if (ExtractFieldsOnly)
  {
    if (Grids[ExtractLevel]->WriteGridX(output, ExtractOutputName, 0) == FAIL) {
      fprintf(stderr, "Error in grid->WriteGridX.\n");
      my_exit(EXIT_FAILURE);
    }
  }
  else
  {
    if (Grids[ExtractLevel]->WriteGrid(output, ExtractOutputName, 0) == FAIL) {
      fprintf(stderr, "Error in grid->WriteGrid.\n");
      my_exit(EXIT_FAILURE);
    }
  }
 
  /* If using particles, open a file to output particle data. */
 
  FILE *StarFile = NULL;
 
  if (StarParticleCreation)
    StarFile = fopen("extract.stars", "wb");
  if (StarFile != NULL) {
    if (fwrite( (void*) &NumberOfStarParticles, sizeof(int), 1, StarFile) != 1)
      perror("error in fwrite1");
    int tempint = NumberOfParticleAttributes + 4;
    fwrite( (void*) &tempint, sizeof(int), 1, StarFile);
    Grids[ExtractLevel]->OutputStarParticleInformation(StarFile);
    fclose(StarFile);
  }
 
  /* --------------------------------------------------------------- */
  /* Write out a grid plot file. */
 
  if ((output = fopen("amr.grid", "w")) == NULL) {
    fprintf(stderr, "Error opening amr.grid for writing.\n");
    my_exit(EXIT_FAILURE);
  }
 
  /* Loop over all levels. */
 
  int Level, grid, Rank, Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  for (Level = 0; Level < MAX_DEPTH_OF_HIERARCHY; Level++) {
 
    /* Loop over this level. */
 
    Temp = LevelArray[Level];
    grid = 0;
    while (Temp != NULL) {
 
      /* Get grid info. */
 
      Temp->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
 
      /* Write out grid data. */
 
      fprintf(output, "Level %"ISYM"   Grid %"ISYM"\n", Level, grid);
 
#define CORNERS
 
#ifdef LINES
      for (i = 0; i < NumberOfLinks[Rank-1]; i++) {
	for (dim = 0; dim < Rank; dim++)
	  if (LinkSide[i][dim] == 0)
	    fprintf(output, "%"FSYM" ", Left[dim]);
	  else
	    fprintf(output, "%"FSYM" ", Right[dim]);
	fprintf(output, "\n");
      }
#endif /* LINES */
 
#ifdef CORNERS
//      for (dim = 0; dim < Rank; dim++)
//	fprintf(output, "%"FSYM" %"FSYM"\n", Left[dim], Right[dim]);
      WriteListOfFloats(output, Rank, Left);
      WriteListOfFloats(output, Rank, Right);
//      fprintf(output, "%"FSYM" %"FSYM" %"FSYM"\n", Right[0]-Left[0], Right[1]-Left[1],
//	      Right[2]-Left[2]);
#endif /* CORNERS */
 
 
      /* Next grid. */
 
      Temp = Temp->NextGridThisLevel;
      grid++;
 
    } // end loop over this level
 
  }
 
  /* Close file. */
 
  fclose(output);
 
  my_exit(EXIT_SUCCESS);
}
