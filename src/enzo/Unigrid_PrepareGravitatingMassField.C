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
      ENZO_FAIL("Error in grid->InitializeGravitatingMassField.\n");
    }
    CurrentGrid->ClearGravitatingMassField();
//  fprintf(stderr, "  PGMF - Initialize & Clear GravitatingMassField\n");
  }
 
  /* Baryons: copy parent density (no interpolation) to regions in
     GravitatingMassField which are beyond the boundary of the current grid. */
 
  if (Grid->ParentGrid != NULL)
   if (CurrentGrid->CopyParentToGravitatingFieldBoundary(
				         Grid->ParentGrid->GridData) == FAIL) {
     ENZO_FAIL("Error in grid->CopyParentToGravitatingFieldBoundary.\n");
   }
 
  /* Baryons: deposit mass into GravitatingMassField. */
 
  //  if (CurrentGrid->AddBaryonsToGravitatingMassField() == FAIL) {
 
//  fprintf(stderr, "  PGMF - DepositBaryons\n");
 
  if (DepositBaryons(Grid, When) == FAIL) {
    ENZO_FAIL("Error in grid->AddBaryonsToGravitatingMassField\n");
  }
 
  /* Particles: go through all the other grids on this level and add all
     their overlapping GravitatingMassFieldParticles to this grid's
     GravitatingMassField.  Handle periodicity properly. */
 
//  fprintf(stderr, "  PGMF - CopyOverlappingParticleMassField\n");
 
  if (CopyOverlappingParticleMassFields(CurrentGrid, MetaData,
					LevelArray, level) == FAIL) {
    ENZO_FAIL("Error in CopyOverlappingParticleMassFields.\n");
  }
 
#ifdef UNUSED
  FLOAT Zero[] = {0,0,0};
  if (CurrentGrid->AddOverlappingParticleMassField(CurrentGrid,Zero) == FAIL) {
    ENZO_FAIL("Error in grid->AddOverlappingParticleMassField.\n");
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
      ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
    }
  }
#endif /* UNUSED */
 
  /* If we are using comoving coordinates, we must adjust the source term. */
 
  if (CommunicationDirection != COMMUNICATION_SEND) {
 
    if (ComovingCoordinates)
      if (CurrentGrid->ComovingGravitySourceTerm() == FAIL) {
	ENZO_FAIL("Error in grid->ComovingGravitySourceTerm.\n");
      }
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  /* Make a guess at the potential for this grid. */
 
  if (level > 0)
    if (CurrentGrid->PreparePotentialField(Grid->ParentGrid->GridData)
	== FAIL) {
      ENZO_FAIL("Error in grid->PreparePotential.\n");

    }
 
 
  return SUCCESS;
}
