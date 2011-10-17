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
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "communication.h"

/* function prototypes */

#ifndef FAST_SIB
int CopyOverlappingParticleMassFields(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
#endif
int DepositBaryons(HierarchyEntry *Grid, FLOAT When);
 
/* EvolveHierarchy function */
 

int PrepareGravitatingMassField1(HierarchyEntry *Grid)
{

  /* declarations */

  int RefinementFactor = RefineBy;
  grid *CurrentGrid = Grid->GridData;

  /* Gravity: initialize and clear the gravitating mass field. */

  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
    if (CurrentGrid->InitializeGravitatingMassField(RefinementFactor) == FAIL){
      ENZO_FAIL("Error in grid->InitializeGravitatingMassField.\n");
    }
    CurrentGrid->ClearGravitatingMassField();
  }

  /* Baryons: copy parent density (no interpolation) to regions in
     GravitatingMassField which are beyond the boundary of the current grid. */

  int CommunicationReceiveIndexLast = CommunicationReceiveIndex;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (Grid->ParentGrid != NULL)
   if (CurrentGrid->CopyParentToGravitatingFieldBoundary(
				         Grid->ParentGrid->GridData) == FAIL) {
     ENZO_FAIL("Error in grid->CopyParentToGravitatingFieldBoundary.\n");
   }
  //  if (CommunicationReceiveIndex != CommunicationReceiveIndexLast)
  //    CommunicationReceiveCurrentDependsOn = CommunicationReceiveIndex-1;

  return SUCCESS;
}

/************************************************************************/

#ifdef FAST_SIB
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, int grid1,
				 SiblingGridList SiblingList[],
				 TopGridData *MetaData, int level,
				 FLOAT When)
#else
int PrepareGravitatingMassField2a(HierarchyEntry *Grid, TopGridData *MetaData,
				 LevelHierarchyEntry *LevelArray[], int level,
				 FLOAT When)
#endif
{
 
  /* declarations */
 
  int grid2;
  grid *CurrentGrid = Grid->GridData;
 
  /* Baryons: deposit mass into GravitatingMassField. */
 
  // IF STATEMENT HERE TO MAKE IT SO NO GAS CONTRIBUTES TO GRAVITY
  if(!SelfGravityGasOff){
    if (DepositBaryons(Grid, When) == FAIL) {
      ENZO_FAIL("Error in DepositBaryons\n");
      printf("      Potential calculated for the gas\n");
    }
  }
 
  /* Particles: go through all the other grids on this level and add all
     their overlapping GravitatingMassFieldParticles to this grid's
     GravitatingMassField.  Handle periodicity properly. */
 
//  fprintf(stderr, "  PGMF - CopyOverlappingParticleMassField\n");
 
#ifdef FAST_SIB
  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
    if (CurrentGrid->CheckForOverlap(SiblingList[grid1].GridList[grid2],
                                     MetaData->LeftFaceBoundaryCondition,
                                     MetaData->RightFaceBoundaryCondition,
                                     &grid::AddOverlappingParticleMassField)
        == FAIL) {
      fprintf(stderr, "Error in grid->AddOverlappingParticleMassFields.\n");
    }
  //  for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
  //fprintf(stderr,"grid %i %i\n", grid1, grid2);

#else
  if (CopyOverlappingParticleMassFields(CurrentGrid, MetaData,
                                        LevelArray, level) == FAIL) {
    ENZO_FAIL("Error in CopyOverlappingParticleMassFields.\n");
  }
#endif

  /* If we are using comoving coordinates, we must adjust the source term. */
 
  if (CommunicationDirection == COMMUNICATION_SEND ||
      CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {

    if (ComovingCoordinates)
      if (CurrentGrid->ComovingGravitySourceTerm() == FAIL) {
	ENZO_FAIL("Error in grid->ComovingGravitySourceTerm.\n");
      }
 
  } // end: if (CommunicationDirection != COMMUNICATION_SEND)
 
  return SUCCESS;
}

/************************************************************************/

int PrepareGravitatingMassField2b(HierarchyEntry *Grid, int level)
{
 
  /* declarations */
 
  grid *CurrentGrid = Grid->GridData;

  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  if (level > 0)

    CurrentGrid->PreparePotentialField(Grid->ParentGrid->GridData);
 
  return SUCCESS;
}
