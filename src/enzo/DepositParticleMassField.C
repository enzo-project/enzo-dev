/***********************************************************************
/
/  DEPOSIT PARTICLES INTO PARTICLEMASSFIELD IN THIS GRID AND ALL SUBGRIDS
/
/  written by: Greg Bryan
/  date:       May, 1995
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
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT Time);
 
int DepositParticleMassField(HierarchyEntry *Grid, FLOAT TimeMidStep)
{
 
  /* Get the time and dt for this grid.  Compute time+1/2 dt. */
 
  if (TimeMidStep < 0)
    TimeMidStep =     Grid->GridData->ReturnTime() +
                  0.5*Grid->GridData->ReturnTimeStep();
 
  /* Initialize the gravitating mass field only if in send-receive mode
     (i.e. this routine is called only once) or if in the first of the
     three communication modes (post-receive). */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
 
    /* Initialize the gravitating mass field parameters (if necessary). */
 
    if (Grid->GridData->InitializeGravitatingMassFieldParticles(RefineBy)
                                                                  == FAIL) {
      ENZO_FAIL("Error in grid->InitializeGravitatingMassFieldParticles.\n");
    }
 
    /* Clear the GravitatingMassFieldParticles. */
 
    if (Grid->GridData->ClearGravitatingMassFieldParticles() == FAIL) {
      ENZO_FAIL("Error in grid->ClearGravitatingMassFieldParticles.\n");
    }
 
//  fprintf(stderr, "--DepositParticleMassField (Send) Initialize & Clear\n");
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  /* Deposit particles to GravitatingMassFieldParticles in this grid. */
 
//  fprintf(stderr, "--DepositParticleMassField Call DepositParticlePositions\n");
 
  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				 GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
 
  /* Recursively deposit particles in children (at TimeMidStep). */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(Grid, Grid->NextGridNextLevel,
					 TimeMidStep)
	== FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren.\n");
    }
 
  return SUCCESS;
}
 
 
 
 
int DepositParticleMassFieldChildren(HierarchyEntry *DepositGrid,
				     HierarchyEntry *Grid, FLOAT DepositTime)
{
 
  /* Deposit particles in Grid into DepositGrid at the given time. */
 
  if (Grid->GridData->DepositParticlePositions(DepositGrid->GridData,
		     DepositTime, GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }
 
  /* Next grid on this level. */
 
  if (Grid->NextGridThisLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridThisLevel,
					 DepositTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(1).\n");
    }
 
  /* Recursively deposit particles in children. */
 
  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(DepositGrid, Grid->NextGridNextLevel,
					 DepositTime) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren(2).\n");

    }
 
 
  return SUCCESS;
}


int DepositParticleMassFieldWithParent(HierarchyEntry *Grid, TopGridData *MetaData,
                                       ChainingMeshStructure *ChainingMesh,
				       FLOAT TimeMidStep)
{

  SiblingGridList SiblingList;
  int grid1;

  /* Get the time and dt for this grid.  Compute time+1/2 dt. */

  if (TimeMidStep < 0)
    TimeMidStep =     Grid->GridData->ReturnTime() +
                  0.5*Grid->GridData->ReturnTimeStep();

  /* Initialize the gravitating mass field only if in send-receive mode
     (i.e. this routine is called only once) or if in the first of the
     three communication modes (post-receive). */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

    /* Initialize the gravitating mass field parameters (if necessary). */

    if (Grid->GridData->InitializeGravitatingMassFieldParticles(RefineBy) == FAIL) {
      ENZO_FAIL("Error in grid->InitializeGravitatingMassFieldParticles.\n");
    }

    /* Clear the GravitatingMassFieldParticles. */

    if (Grid->GridData->ClearGravitatingMassFieldParticles() == FAIL) {
      ENZO_FAIL("Error in grid->ClearGravitatingMassFieldParticles.\n");
    }
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)

  /* Deposit particles to GravitatingMassFieldParticles in this grid. */

  if (Grid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
                                 GRAVITATING_MASS_FIELD_PARTICLES) == FAIL) {
    ENZO_FAIL("Error in grid->DepositParticlePositions.\n");
  }

  /* Recursively deposit particles in children (at TimeMidStep). */

  if (Grid->NextGridNextLevel != NULL)
    if (DepositParticleMassFieldChildren(Grid, Grid->NextGridNextLevel,
                                         TimeMidStep) == FAIL) {
      ENZO_FAIL("Error in DepositParticleMassFieldChildren.\n");
    }

  /* Copy particles from parent and siblings from the upper level */

  if (Grid->ParentGrid != NULL) {

    /* First do the parent grid */

    if (Grid->ParentGrid->GridData->DepositParticlePositions(Grid->GridData, TimeMidStep,
				     GRAVITATING_MASS_FIELD_PARTICLES) == FAIL)
      ENZO_FAIL("Error in grid->DepositParticlePositions: ParentGrid.\n");

    /* Now do the siblings
       Based on CommunicationTransferSubgridParticles.C
    */

    // Get a list of possible siblings from the chaining mesh

    Grid->GridData->FastSiblingLocatorFindSiblings
      (ChainingMesh, &SiblingList, MetaData->LeftFaceBoundaryCondition,
       MetaData->RightFaceBoundaryCondition);

    // Deposit particles from sibling

    for (grid1 = 0; grid1 < SiblingList.NumberOfSiblings; grid1++)
      if (Grid->ParentGrid->GridData != SiblingList.GridList[grid1])
	if (SiblingList.GridList[grid1]->DepositParticlePositions(
		  Grid->GridData, TimeMidStep, GRAVITATING_MASS_FIELD_PARTICLES
                                                                    ) == FAIL)
	  ENZO_FAIL("Error in grid->DepositParticlePositions: Sibling.\n");

    // Cleanup

    delete [] SiblingList.GridList;

  }

  return SUCCESS;
}


