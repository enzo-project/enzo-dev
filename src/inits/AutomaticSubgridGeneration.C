/***********************************************************************
/
/  AUTOMATICALLY GENERATES A SET OF REFINED INITIAL GRIDS
/
/  written by: Greg Bryan
/  date:       June, 2005
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Parameters.h"

/* function prototypes */

int GenerateRealization(parmstruct *Parameters, parmstruct *SubGridParameters);
int intpow(int base, int exponent);



int AutomaticSubgridGeneration(parmstruct *Parameters)
{

  /* locals */

  int level, dim, StartIndex[3], EndIndex[3], LevelDims[3];
  parmstruct LevelParms[MAX_LEVELS];
  FLOAT FinalLeftEdge[3], FinalRightEdge[3];

  /* Some error checking. */

  if (Parameters->StartIndexInNewCenterTopGridSystem[0] != INT_UNDEFINED) {
    fprintf(stderr, "Cannot set both StartIndexInNewCenterTopGridSystem and MaximumInitialRefinementLevel; pick only one.\n");
    exit(EXIT_FAILURE);
  }

  if (Parameters->MaximumInitialRefinementLevel >= MAX_LEVELS) {
    fprintf(stderr, "You must increase the value of MAX_LEVELS\n");
    exit(EXIT_FAILURE);
  }

  /* Loop over levels to be generated. */

  for (level=0; level <= Parameters->MaximumInitialRefinementLevel; level++) {

    /* create a copy of the parameter file. */

    LevelParms[level] = *Parameters;

    /* compute GridDims and StartIndex (except for level 0). */

    if (level > 0) {
      for (dim = 0; dim < Parameters->Rank; dim++) {

	/* Compute the dimensions at the level which is one coarser so that we 
	   start and end of a grid boundary of that level's grid. */

	LevelDims[dim] = Parameters->RootGridDims[dim]*
	                 intpow(Parameters->RefineBy, level-1);
	
	/* Compute the start and end indices in the coarse level grid.   This is the
	   bit of code which actually sets the size of the regions. */

	StartIndex[dim] = max(int(LevelDims[dim]*Parameters->RefineRegionLeftEdge[dim]) - Parameters->AutomaticSubgridBuffer, 0);
	EndIndex[dim] = int(LevelDims[dim]*Parameters->RefineRegionRightEdge[dim]) + Parameters->AutomaticSubgridBuffer;

	/* Convert to the most refined level for start index. */

	LevelParms[level].StartIndex[dim] = StartIndex[dim]*intpow(Parameters->RefineBy, Parameters->MaximumInitialRefinementLevel-(level-1));

	/* Compute the dimensions of the grid on this level. */

	if (Parameters->GridRefinement)
	  LevelParms[level].GridDims[dim] = 
	    (EndIndex[dim] - StartIndex[dim] + 1)*Parameters->RefineBy;
	if (Parameters->ParticleRefinement)
	  LevelParms[level].ParticleDims[dim] = 
	    (EndIndex[dim] - StartIndex[dim] + 1)*Parameters->RefineBy;

      }	 // end: loop over dims
    } // end: if (level > 0)

    /* Set the names. */

    if (Parameters->ParticlePositionName != NULL)
      sprintf(LevelParms[level].ParticlePositionName = new char[MAX_LINE_LENGTH],
	      "%s.%1d", Parameters->ParticlePositionName, level);
    if (Parameters->ParticleVelocityName != NULL)
      sprintf(LevelParms[level].ParticleVelocityName = new char[MAX_LINE_LENGTH],
	      "%s.%1d", Parameters->ParticleVelocityName, level);
    if (Parameters->ParticleMassName != NULL)
      sprintf(LevelParms[level].ParticleMassName = new char[MAX_LINE_LENGTH],
	      "%s.%1d", Parameters->ParticleMassName, level);
    if (Parameters->GridDensityName != NULL)
      sprintf(LevelParms[level].GridDensityName = new char[MAX_LINE_LENGTH],
	      "%s.%1d", Parameters->GridDensityName, level);
    if (Parameters->GridVelocityName != NULL)
      sprintf(LevelParms[level].GridVelocityName = new char[MAX_LINE_LENGTH],
	      "%s.%1d", Parameters->GridVelocityName, level);

    /* Set Grid and Particle Refinement. */

    for (dim = 0; dim < Parameters->Rank; dim++) {
      LevelParms[level].GridRefinement     = Parameters->GridRefinement*
           intpow(Parameters->RefineBy, Parameters->MaximumInitialRefinementLevel-level);
      LevelParms[level].ParticleRefinement = Parameters->ParticleRefinement*
           intpow(Parameters->RefineBy, Parameters->MaximumInitialRefinementLevel-level);
    } // end: loop over dims

    printf("level %"ISYM"   ", level);
    if (Parameters->GridRefinement) {
      printf("GridDims = ");
      for (dim = 0; dim < Parameters->Rank; dim++)
	printf("%"ISYM" ", LevelParms[level].GridDims[dim]);
    }
    if (Parameters->ParticleRefinement) {
      printf("   ParticleDims = ");
      for (dim = 0; dim < Parameters->Rank; dim++)
	printf("%"ISYM" ", LevelParms[level].ParticleDims[dim]);
    }
    printf("\n");

  } // end: loop over levels

  /* Create the output for the enzo parameter file. */

  FILE *fptr = fopen("EnzoMultigridParameters", "w");
  if (fptr != NULL) {
    fprintf(fptr, "CosmologySimulationNumberOfInitialGrids = %"ISYM"\n",
	    Parameters->MaximumInitialRefinementLevel+1);
    for (level=1; level <= Parameters->MaximumInitialRefinementLevel; level++) {

      /* Compute the position of the left corner of the grid in the final, moved
	 coordinate system.   This is only required for output here and should
	 match the code in GenerateRealization. */

      for (dim = 0; dim < Parameters->Rank; dim++) {
	if (LevelParms[level].NewCenter[0] != INT_UNDEFINED && 
	    LevelParms[level].MaxDims[dim] != LevelParms[level].ParticleDims[dim]*
	                              LevelParms[level].ParticleRefinement)
	  FinalLeftEdge[dim] = 
	    FLOAT((LevelParms[level].StartIndex[dim] + LevelParms[level].MaxDims[dim]/2-1 -
		 LevelParms[level].NewCenter[dim] + LevelParms[level].MaxDims[dim]      ) 
		% LevelParms[level].MaxDims[dim]) /
	    FLOAT(LevelParms[level].MaxDims[dim]);
	else
	  FinalLeftEdge[dim] = FLOAT(LevelParms[level].StartIndex[dim]) / 
	                       FLOAT(LevelParms[level].MaxDims[dim]);
	FinalRightEdge[dim] = FinalLeftEdge[dim] + 
	  FLOAT(LevelParms[level].ParticleDims[dim]*LevelParms[level].ParticleRefinement)/
	  FLOAT(LevelParms[level].MaxDims[dim]);
      } // end: loop over dims

      /* Output parameters to a fragment of an enzo parameter file. */

      fprintf(fptr, "CosmologySimulationGridDimension[%"ISYM"]     = ", level);
      for (dim = 0; dim < Parameters->Rank; dim++)
	if (Parameters->GridRefinement)
	  fprintf(fptr, "%"ISYM" ", LevelParms[level].GridDims[dim]);
	else
	  fprintf(fptr, "%"ISYM" ", LevelParms[level].ParticleDims[dim]);
      fprintf(fptr, "\n");

      fprintf(fptr, "CosmologySimulationGridLeftEdge[%"ISYM"]     = ", level);
      for (dim = 0; dim < Parameters->Rank; dim++)
	fprintf(fptr, "%"FSYM" ", FinalLeftEdge[dim]);
      fprintf(fptr, "\n");

      fprintf(fptr, "CosmologySimulationGridRightEdge[%"ISYM"]     = ", level);
      for (dim = 0; dim < Parameters->Rank; dim++)
	fprintf(fptr, "%"FSYM" ", FinalRightEdge[dim]);
      fprintf(fptr, "\n");

      fprintf(fptr, "CosmologySimulationGridLevel[%"ISYM"]         = %"ISYM"\n", level, level);

    } // end: loop over levels

    fclose(fptr);

  } // end: if (fptr != NULL)

  /* Now, actually generate each level. */

  for (level = 0; level < Parameters->MaximumInitialRefinementLevel; level++) {

    printf("Generating level %"ISYM"\n", level);
    if (GenerateRealization(LevelParms+level, LevelParms+level+1) 
	== FAIL) {
      fprintf(stderr, "Error in GenerateRealization on level %"ISYM"\n", level);
      exit(EXIT_FAILURE);
    }

  } // end: loop over levels

  printf("Generating level %"ISYM"\n", level);
  GenerateRealization(LevelParms+Parameters->MaximumInitialRefinementLevel, NULL);

  return SUCCESS;
}

/* Integer power, this isn't efficient but it keeps things integer. */

int intpow(int base, int exponent)
{
  int ret = 1;
  for (int i = 0; i < exponent; i++)
    ret *= base;

  return ret;
}
