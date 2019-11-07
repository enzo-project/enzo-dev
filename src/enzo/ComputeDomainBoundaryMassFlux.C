/**********************************************************************
/
/ COMPUTE DOMAIN BOUNDARY MASS FLUX
/
/ written by: Andrew Emerick
/ date:       January 2017
/ modified1:
/
/ Purpose:
/   This routine sums up the external boundary mass flux for each
/   field assigned in Grid_PrepareBoundaryMassFluxFieldNumbers.
/
/   NOTE: If no outflow boundaries are set, this will just give
/         zeros for all fields. Grid_ComputeDomainBoundaryMassFlux assumes
/         mass flux for each field, reporting the mass that leaves
/         the grid in solar masses. However, this could (in principle)
/         be generalized to compute flux of any quantity in any units
/
/*********************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"


int ComputeDomainBoundaryMassFlux(HierarchyEntry *Grids[], int level,
                                  int NumberOfGrids,
                                  TopGridData *MetaData){

  if (!StoreDomainBoundaryMassFlux) return SUCCESS;

  if (level != 0) return SUCCESS; // only do for root grid

  int grid1;

  float allgrid_BoundaryMassFluxContainer[MAX_NUMBER_OF_BARYON_FIELDS];
  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++) allgrid_BoundaryMassFluxContainer[i] = 0.0;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++){

    if( Grids[grid1]->GridData->ComputeDomainBoundaryMassFlux(allgrid_BoundaryMassFluxContainer,
                                                              MetaData) == FAIL)
      ENZO_FAIL("Error in grid->ComputeDomainBoundaryMassFlux.\n");
  }

  /* now communicate */
  CommunicationSumValues(allgrid_BoundaryMassFluxContainer, MAX_NUMBER_OF_BARYON_FIELDS);

  /* now store total in root grid for output */
  if (MyProcessorNumber == ROOT_PROCESSOR){

    FILE *fptr;

    /* check if this is the first time through and write header */
    if (BoundaryMassFluxContainer[0] <= 0){
      if ((fptr = fopen(BoundaryMassFluxFilename, "w")) == NULL){
        ENZO_VFAIL("Error opening file %s\n", BoundaryMassFluxFilename); 
      }
      fprintf(fptr, "Cycle_Number Time ");
      for(int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++){
        if(BoundaryMassFluxFieldNumbers[i] < 0) break;
        fprintf(fptr, " %s ", DataLabel[ BoundaryMassFluxFieldNumbers[i] ]);
      }
      fprintf(fptr, "\n");
      fclose(fptr);
    }

    /* now add to the cumulative boundary mass flux container */
    for( int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++){
        BoundaryMassFluxContainer[i] += max(allgrid_BoundaryMassFluxContainer[i], 0.0);
    }

    /* now output the amount of mass outflow in this cycle */
    if ((fptr = fopen(BoundaryMassFluxFilename, "a")) == NULL){
      ENZO_VFAIL("Error opening file %s\n", BoundaryMassFluxFilename);
    }

    fprintf(fptr, "%"ISYM" %"GSYM" ", MetaData->CycleNumber, MetaData->Time);
    for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++){
      if( BoundaryMassFluxFieldNumbers[i] < 0) break;
      fprintf(fptr, "%"GSYM" ", allgrid_BoundaryMassFluxContainer[i]);
    }
    fprintf(fptr, "\n");
    fclose(fptr);

  }

  return SUCCESS;
}
