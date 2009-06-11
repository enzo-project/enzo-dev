/***********************************************************************
/
/  PREPARE THE GRAVITATING MASS FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
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
 
int CopyOverlappingParticleMassFields(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
int DepositBaryons(HierarchyEntry *Grid, FLOAT When);
 
/* EvolveHierarchy function */
 
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level,
                                FLOAT When)
{
 
  /* declarations */
 
  int RefinementFactor = RefineBy;
  grid *CurrentGrid = Grid->GridData;
 
  /* Gravity: initialize and clear the gravitating mass field. */
 
  if (CommunicationDirection != COMMUNICATION_RECEIVE) {
    if (CurrentGrid->InitializeGravitatingMassField(RefinementFactor) == FAIL){
      fprintf(stderr, "Error in grid->InitializeGravitatingMassField.\n");
      ENZO_FAIL("");
    }
    CurrentGrid->ClearGravitatingMassField();
//  fprintf(stderr, "  PGMF - Initialize & Clear GravitatingMassField\n");
  }
 
  /* Baryons: copy parent density (no interpolation) to regions in
     GravitatingMassField which are beyond the boundary of the current grid. */
 
  if (Grid->ParentGrid != NULL)
   if (CurrentGrid->CopyParentToGravitatingFieldBoundary(
				         Grid->ParentGrid->GridData) == FAIL) {
     fprintf(stderr, "Error in grid->CopyParentToGravitatingFieldBoundary.\n");
     ENZO_FAIL("");
   }
 
  /* Baryons: deposit mass into GravitatingMassField. */
 
  //  if (CurrentGrid->AddBaryonsToGravitatingMassField() == FAIL) {
 
//  fprintf(stderr, "  PGMF - DepositBaryons\n");
 
  if (DepositBaryons(Grid, When) == FAIL) {
    fprintf(stderr, "Error in grid->AddBaryonsToGravitatingMassField\n");
    ENZO_FAIL("");
  }
 
  /* Particles: go through all the other grids on this level and add all
     their overlapping GravitatingMassFieldParticles to this grid's
     GravitatingMassField.  Handle periodicity properly. */
 
//  fprintf(stderr, "  PGMF - CopyOverlappingParticleMassField\n");
 
  if (CopyOverlappingParticleMassFields(CurrentGrid, MetaData,
					LevelArray, level) == FAIL) {
    fprintf(stderr, "Error in CopyOverlappingParticleMassFields.\n");
    ENZO_FAIL("");
  }
 
#ifdef UNUSED
  FLOAT Zero[] = {0,0,0};
  if (CurrentGrid->AddOverlappingParticleMassField(CurrentGrid,Zero) == FAIL) {
    fprintf(stderr, "Error in grid->AddOverlappingParticleMassField.\n");
    ENZO_FAIL("");
  }
#endif /* UNUSED */
 
  /* Particles: deposit particles in the parent grid into GravitatingMassField
     (they should only be in the boundaries).  We should really do this for
     all parent grids.  Or better yet do a PP summation. sigh. */
 
#ifdef UNUSED
  if (Grid->ParentGrid != NULL) {
 
/* The following is done to allow this to be parallelized. */
//    FLOAT TimeMidStep = CurrentGrid->ReturnTime() +
//                        When*CurrentGrid->ReturnTimeStep();
      FLOAT TimeMidStep = Grid->ParentGrid->GridData->ReturnTime();
 
    if (Grid->ParentGrid->GridData->DepositParticlePositions(CurrentGrid,
			       TimeMidStep, GRAVITATING_MASS_FIELD) == FAIL) {
      fprintf(stderr, "Error in grid->DepositParticlePositions.\n");
      ENZO_FAIL("");
    }
  }
#endif /* UNUSED */
 
  /* If we are using comoving coordinates, we must adjust the source term. */
 
  if (CommunicationDirection != COMMUNICATION_SEND) {
 
    if (ComovingCoordinates)
      if (CurrentGrid->ComovingGravitySourceTerm() == FAIL) {
	fprintf(stderr, "Error in grid->ComovingGravitySourceTerm.\n");
	ENZO_FAIL("");
      }
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  /* Make a guess at the potential for this grid. */
 
  if (level > 0)
    if (CurrentGrid->PreparePotentialField(Grid->ParentGrid->GridData)
	== FAIL) {
      fprintf(stderr, "Error in grid->PreparePotential.\n");
      ENZO_FAIL("");
    }
 
 
  return SUCCESS;
}
