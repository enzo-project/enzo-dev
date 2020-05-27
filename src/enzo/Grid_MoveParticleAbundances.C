
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

int DetermineNumberOfAbundanceAttributes(void);

int grid::MoveParticleAbundances(int NumberOfGrids, grid* FromGrid[]){

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }

  /* Determine total number of local particles. */

  int NumberOfSubgridParticles = 0;
  int TotalNumberOfParticles = NumberOfParticles;
  int i, j, grid;

  for (grid = 0; grid < NumberOfGrids; grid++)
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
      NumberOfSubgridParticles += FromGrid[grid]->NumberOfParticles;
  if (NumberOfSubgridParticles == 0)
    return SUCCESS;

  TotalNumberOfParticles += NumberOfSubgridParticles;

  this->AllocateStellarAbundances(NumberOfSubgridParticles);

  int Index = 0;

  int num_abundances = DetermineNumberOfAbundanceAttributes();

  for (grid = 0; grid < NumberOfGrids; grid++){
    for (j = 0; j < num_abundances; j++){
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
