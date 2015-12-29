/***********************************************************************
/
/  GRID CLASS (FIND AND FLAG PARTICLES AS MUST REFINE FROM LIST OF PARTICLE IDS)
/
/  written by: Christine Simpson & Greg Bryan
/  date:       December 2009
/  modified1:
/
/  PURPOSE:Find Particles which match the provided list of particle id 
/  numbers and flag them as must refine particles.  The particle list must 
/  be called MustRefineParticlesFlaggingList.in.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int grid::MustRefineParticlesFlagFromList()
{
  
  int i, dim, LowerIndex, UpperIndex, MidPoint, ParticlesFound, ParticlesFlagged, *ParticleNumberList;
  int NumberOfParticles;
  FILE *fptr;
  float tempf;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /*Read in particle indicies*/

  fptr=fopen("MustRefineParticlesFlaggingList.in","r");
 
  char line[MAX_LINE_LENGTH];
  ParticlesFound = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#') ParticlesFound++;
   
  if (debug) printf("P(%d): ParticlesFound = %d\n",MyProcessorNumber,ParticlesFound);
  
  ParticleNumberList=new int[ParticlesFound]; //to allocate array space

  rewind(fptr);
  i = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#')
      if (sscanf(line, "%"ISYM, &ParticleNumberList[i]) == 1)
        i++;

  fclose(fptr);

  /* Loop over particles. A bisection search is employed to speed things up.
     For this to work, the ParticleNumberList must be sorted before it is 
     read in. */
  ParticlesFlagged=0.;
  for (i = 0; i < NumberOfParticles; i++) {
    /* Find place in sorted ParticleNumberList by bisection. */

    LowerIndex = -1;
    UpperIndex = ParticlesFound;

    while (UpperIndex-LowerIndex > 1) {
      MidPoint = (LowerIndex + UpperIndex)/2;
      if (ParticleNumber[i] > ParticleNumberList[MidPoint])
	LowerIndex = MidPoint;
      else
	UpperIndex = MidPoint;
    }
    
    /* If found, set left/right edge and output */
    if (ParticleNumber[i] == ParticleNumberList[LowerIndex] ||
	ParticleNumber[i] == ParticleNumberList[UpperIndex]  ) {
      ParticleType[i] = PARTICLE_TYPE_MUST_REFINE;
      ParticlesFlagged=ParticlesFlagged+1.;
     } // end: if (index != -1)

  } // end: loop over particles
  if (ParticlesFlagged > 0)
    printf("P(%d): MustRefineParticlesFlagFromList:ParticlesFlagged=%d\n",MyProcessorNumber,ParticlesFlagged);
  delete [] ParticleNumberList;
  return SUCCESS;
}
