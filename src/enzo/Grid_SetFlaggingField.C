/***********************************************************************
/
/  GRID CLASS (CLEAR THE FLAGGING FIELD)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, Aug. 2004: added refinement by shear.
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* The following is defined in Grid_DepositParticlePositions.C. */
 
extern float DepositParticleMaximumParticleMass;
 
 
int grid::SetFlaggingField(int &NumberOfFlaggedCells, int level)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  NumberOfFlaggedCells = INT_UNDEFINED;
 
  /***********************************************************************/
  /* beginning of Cell flagging criterion routine                        */

  int method;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
 
  switch (CellFlaggingMethod[method]) {
 
  case 0:   /* no action */
    NumberOfFlaggedCells = (NumberOfFlaggedCells == INT_UNDEFINED ?
			    0 : NumberOfFlaggedCells);
    break;
 
    /* ==== METHOD 1: BY SLOPE ==== */
 
  case 1:
 
    /* flag all points needing extra resolution (FlagCellsToBeRefinedBySlop
       returns the number of flagged cells). */
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedBySlope();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedBySlope.");
    }
    break;
 
    /* ==== METHOD 2: BY BARYON MASS OR OVERDENSITY ==== */
 
  case 2:
 
    /* allocate and clear mass flagging field */
 
    this->ClearMassFlaggingField();
 
    /* baryons: add baryon density to mass flagging field (so the mass
       flagging field contains the mass in the cell (not the density) */
 
    if (this->AddFieldMassToMassFlaggingField() == FAIL) {
      ENZO_FAIL("Error in grid->AddFieldMassToMassFlaggingField.");
    }
 
    /* flag all points that need extra resolution (FlagCellsToBeRefinedByMass
       return the number of flagged cells). */
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByMass (2).");
    }
    break;
 
    /* ==== METHOD 3: BY SHOCKS ==== */
 
  case 3:
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShocks();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByShocks.");
    }
    break;
 
    /* ==== METHOD 4: BY PARTICLE MASS ==== */
 
  case 4:

    /* All of the calculation of particle mass flagging fields are
       done in grid::SetParticleMassFlaggingField now. */
 
    /* Flag all points that need extra resolution (FlagCellsToBeRefinedByMass
       return the number of flagged cells). */
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMass(level, method);
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByMass (4).");
    }
    break;
 
    /* ==== METHOD 5: (disabled)  ==== */

    /* ==== METHOD 6: BY JEANS LENGTH ==== */
 
  case 6:
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByJeansLength();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByJeansLength.");
    }
    break;
 
    /* ==== METHOD 7: BY COOLING TIME < DX/SOUND SPEED ==== */
 
  case 7:
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByCoolingTime();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByCoolingTime.");
    }
    break;
 
    /* ==== METHOD 8: BY POSITION OF MUST-REFINE PARTICLES  ==== */
 
  case 8:

    /* Searching for must-refine particles now done in
       grid::SetParticleMassFlaggingField and stored in
       ParticleMassFlaggingField.  This is checked in method #4, which
       is automatically turned if method #8 is specified. */

    break;
 
    /* ==== METHOD 9: BY SHEAR ==== */
 
  case 9:

    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShear();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByShear.");
    }
    break;

    /* ==== METHOD 10: BY OPTICAL DEPTH ==== */
 
  case 10:
#ifdef TRANSFER
    if (RadiativeTransfer) {
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByOpticalDepth();
      if (NumberOfFlaggedCells < 0) {
	ENZO_FAIL("Error in grid->FlagCellsByOpticalDepth.");
      }
    }
#endif /* TRANSFER */
    break;

    /* ==== METHOD 11: BY RESISTIVE LENGTH ==== */

  case 11:
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByResistiveLength();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByResistiveLength.");
    }
    break;

    /* ==== METHOD 12: FORCE REFINEMENT TO SOME LEVEL IN A SET REGION ==== */
 
  case 12:
    if (level < MustRefineRegionMinRefinementLevel) {
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMustRefineRegion(level);
      if (NumberOfFlaggedCells < 0) {
	ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByMustRefineRegion.");
      }
    }
    break;
 
    /* ==== METHOD 13: FORCE REFINEMENT BASED ON METALLICITY OF GAS ==== */
    
  case 13:
    if (level < MetallicityRefinementMinLevel) {
      NumberOfFlaggedCells = this->FlagCellsToBeRefinedByMetallicity(level);
      if (NumberOfFlaggedCells < 0) {
	ENZO_FAIL("Error in grid->FlagCellsToBeRefinedByMetallicity.");
      }
    }
    break;
    
    /* ==== METHOD 14: Refine around Shockwaves ==== */
  case 14:
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByShockwaves(level);
    if (NumberOfFlaggedCells < 0) {
      fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByShockwaves.\n");
      return FAIL;
    }
    break;

      /* ==== METHOD 15: Refine by Second Derivative ==== */
  case 15:
 
    /* flag all points needing extra resolution 
     * (FlagCellsToBeRefinedBySecondDerivative)
       returns the number of flagged cells). */
 
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedBySecondDerivative();
    if (NumberOfFlaggedCells < 0) {
      ENZO_FAIL("Error in grid->FlagCellsToBeRefinedBySecondDerivative.");
    }
    break;
    
    /* ==== METHOD 16: Refine on total Jeans length ==== */
  case 16: 
    //    this->SolveForPotential(level);
    NumberOfFlaggedCells = this->FlagCellsToBeRefinedByTotalJeansLength();
    if (NumberOfFlaggedCells < 0) {
      fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByTotalJeansLength.\n");
      return FAIL;
    }
    break;

    /* ==== METHOD 17: undefined ==== */

    /* ==== METHOD 18: BY POSITION OF MUST-REFINE PARTICLES ONLY ABOVE A CERTAIN MASS  ==== */
  case 18:

    /* Searching for must-refine particles now done in
       grid::SetParticleMassFlaggingField and stored in
       ParticleMassFlaggingField.  This is checked in method #4, which
       is automatically turned if method #8 is specified. */

    break;
 
    /* ==== METHOD 100: UNDO REFINEMENT IN SOME REGIONS ==== */
 
    /* Must be done last ... */
  case 100:
    this->FlagCellsToAvoidRefinement();
    if (NumberOfFlaggedCells < 0)
      ENZO_FAIL("Error in grid->FlagCellsToAvoidRefinement");
    break;

  case 101:
    this->FlagCellsToAvoidRefinementRegion(level);
    if (NumberOfFlaggedCells < 0)
      ENZO_FAIL("Error in grid->FlagCellsToAvoidRefinementRegion");
    break;

  case INT_UNDEFINED:
    break;
 
  default:
    ENZO_VFAIL("CellFlaggingMethod[%"ISYM"] = %"ISYM" unknown\n", method,
	    CellFlaggingMethod[method])
 
  }

  } // ENDFOR methods
 
  /* End of Cell flagging criterion routine                              */
  /***********************************************************************/
 
  if (NumberOfFlaggedCells == INT_UNDEFINED) {
    ENZO_FAIL("No valid CellFlaggingMethod specified.");
  }
 
#ifdef MPI_INSTRUMENTATION
  counter[4]++;
  timer[4] += NumberOfFlaggedCells;
#endif /* MPI_INSTRUMENTATION */
 
  if (debug1)

    printf("SetFlaggingField[method = %"ISYM"]: NumberOfFlaggedCells = %"ISYM".\n",
	   method, NumberOfFlaggedCells);
 
  return SUCCESS;
 
}
