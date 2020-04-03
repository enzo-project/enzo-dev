 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#ifdef INDIVIDUALSTAR
int grid::MoveParticleAbundances(int NumberOfGrids, grid* FromGrid[]){

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }
 
  /* Determine total number of local particles. */

  int NumberOfSubgridParticles = 0;
  int TotalNumberOfParticles = NumberOfParticles;
  int i, j, grid, dim, *Type;
  PINT *Number;
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
      NumberOfSubgridParticles += FromGrid[grid]->NumberOfParticles;
  if (NumberOfSubgridParticles == 0) 
    return SUCCESS;
  
  TotalNumberOfParticles += NumberOfSubgridParticles;

  this->AllocateStellarAbundances(NumberOfSubgridParticles);

  int Index = 0;
  int offset = 0;

  if (IndividualStarTrackAGBMetalDensity) offset++;  // AGB Metal tracer - ExtraField0
  if (IndividualStarPopIIIFormation)      offset = offset + 2;  // PopIII Metal Tracer - ExtraField1
  if (IndividualStarTrackSNMetalDensity)  offset = offset + 2; // MetalSNIaDensity, MetalSNIIDensity

  for (grid = 0; grid < NumberOfGrids; grid++){
    for (j = 0; j < StellarYieldsNumberOfSpecies + offset; j++){
      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++){
        StellarAbundances[j][Index+i] = FromGrid[grid]->StellarAbundances[j][i];
      }
    }

    Index += FromGrid[grid]->NumberOfParticles;
  }

  /* Do not delete particle information -- */

  return SUCCESS;

}
#endif
