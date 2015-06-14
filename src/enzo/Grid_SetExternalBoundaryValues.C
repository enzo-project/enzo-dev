/***********************************************************************
/
/  GRID CLASS (SET EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Copy the current baryon fields to the old baryon fields/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::SetExternalBoundaryValues(ExternalBoundary *Exterior)
{
  int dim, field;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* For Wave Pool problem, compute the new inflow boundary conditions. */
 
  if (ProblemType == 2)
    if (Exterior->SetWavePoolBoundary(Time) == FAIL) {
      ENZO_FAIL("Error in exterior->SetWavePoolBoundary.\n");
    }
 
  /* For Shock Pool problem, compute the new inflow boundary conditions. */
 
  if (ProblemType == 3)
    if (Exterior->SetShockPoolBoundary(Time) == FAIL) {
      ENZO_FAIL("Error in exterior->SetShockPoolBoundary.\n");
    }

  /* For Galaxy Sim w/ ICM Wind, compute the new inflow boundary conditions. */

  if (ProblemType == 31)
    if (Exterior->SetGalaxySimulationBoundary(Time) == FAIL) {
      ENZO_FAIL("Error in exterior->SetGalaxySimulationBoundary.\n");
    }
 
  /* For the DoubleMach problem, set the bew inflow boundary conditions. */
 
  if (ProblemType == 4)
    if (Exterior->SetDoubleMachBoundary(Time, CellLeftEdge[0], CellWidth[0])
	== FAIL) {
      ENZO_FAIL("Error in exterior->SetDoubleMachBoundary.\n");
    }
 
  /* For 2D/3D Noh problem apply time-dependent boundary conditions on Right faces
     before applying reflecting BCs on the Left ones, thus, taking care of the
     (0,1) and (1,0) corners. */

  if (ProblemType == 9)
    if (this->ComputeExternalNohBoundary() == FAIL) {
      ENZO_FAIL("Error in grid->ComputeExternalNohBoundary.\n");
    }

  /* For the SetWengenCollidingFlowBoundary problem, set the inflow boundary conditions. */
 
  if (ProblemType == 201 && (EOSSoundSpeed > 0)) 
    if (Exterior->SetWengenCollidingFlowBoundary(Time, CellLeftEdge[0], CellWidth[0])
	== FAIL) {
      ENZO_FAIL("Error in exterior->SetWengenCollidingFlowBoundary.\n");
    }


  /* Compute offset from corner of domain. */
 
  int GridOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (dim < GridRank)
      GridOffset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			     CellWidth[dim][0]);
    else
      GridOffset[dim] = 0;
 
  /* loop through fields, setting each */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
 
    if (BaryonField[field] == NULL) {
      ENZO_VFAIL("Baryon field missing. %i %i\n", field, NumberOfBaryonFields)
    }

#ifdef OOC_BOUNDARY
    ExternalBoundaryField = field;
#endif

    if (Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
				      GridStartIndex, GridEndIndex,
				      BaryonField[field], FieldType[field])
	== FAIL) {
      ENZO_FAIL("Error in Exterior->SetExternalBoundary.\n");
    }
 
  }

  /*If there's a magnetic field, set it as well.
    It is still unclear if this is a valid way to do things: 
    The Pseudo-Vector nature of B poses a problem that I haven't sorted out all the way.
    Currently, it's being set as if it were a Plain vector.*/

  if(UseMHDCT && NumberOfBaryonFields > 0)
    {

      if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
				       GridStartIndex, GridEndIndex,
				       MagneticField[0],Bfield1) == FAIL)
	ENZO_FAIL("Error in Exterior->SetMagneticBoundary, B1.");

      if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
				       GridStartIndex, GridEndIndex,
				       MagneticField[1],Bfield2) == FAIL)
	ENZO_FAIL("Error in Exterior->SetMagneticBoundary, B2.");

      if(Exterior->SetMagneticBoundary(GridRank, GridDimension, GridOffset,
				       GridStartIndex, GridEndIndex,
				       MagneticField[2],Bfield3) == FAIL)
	ENZO_FAIL("Error in Exterior->SetMagneticBoundary, B3.");

      //centeredB is set using FieldType == VelocityX since the FieldType array is 
      //compared against.  This will be irrelevant when CenteredB gets stored in BaryonField.
      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
				CenteredB[0], Density) == FAIL)
ENZO_FAIL("Error: Something's wrong with the CenteredB[0] boundary.");
      
      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
				CenteredB[1], Density) == FAIL)
ENZO_FAIL("Error: Something's wrong with the CenteredB[1] boundary.");

      if( Exterior->SetExternalBoundary(GridRank, GridDimension, GridOffset,
					GridStartIndex, GridEndIndex,
				CenteredB[2], Density) == FAIL)
ENZO_FAIL("Error: Something's wrong with the CenteredB[2] boundary.");
 
 
    }// if(UseMHDCT)

  /* Now we handle the particles (if any). */
 
  if (NumberOfParticles > 0)
 
    if (Exterior->SetExternalBoundaryParticles(GridRank, NumberOfParticles,
					       ParticlePosition,
					       ParticleVelocity) == FAIL) {
      ENZO_FAIL("Error in Exterior->SetExternalBoundaryParticles.\n");

    }
 
  return SUCCESS;
 
}
