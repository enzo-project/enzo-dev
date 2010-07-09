#ifndef FLUX_FIX
/***********************************************************************
/
/  GRID CLASS (CORRECT SOLUTION GIVEN ORIGINAL AND REFINED FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       January, 2003
/              Include extra fields beyond Metallicity!
/  modified2: David Collins, 2005
/              Updated algebra so Cosmological Expansion is also
/              conservative.  This fix also came with fixes to euler.src and
/              Grid_GetProjectedBoundaryFluxes.C, so make sure you get those.
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// Given both the original and refined fluxes for a subgrid in the current
//   grid, correct the solution in the cells just outside the solution to
//   reflect the new fluxes (which are determined from the solution of
//   the subgrid).
//   Note that if the subgrid is on the boundary of the current grid, we
//     do not correct the values but instead replace the boundary fluxes
//     for the current time step (BoundaryFluxesThisTimeStep).
 
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::CorrectForRefinedFluxes(fluxes *InitialFluxes,
				  fluxes *RefinedFluxes,
				  fluxes *BoundaryFluxesThisTimeStep)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro)
    return SUCCESS;
 
  /* declarations */
 
  int i1, i2, i, j, k, dim, field, ffield, index;
  int FieldIndex, FluxIndex, GridFluxIndex;
  int CorrectLeftBoundaryFlux, CorrectRightBoundaryFlux, Offset;
  int Dim[MAX_DIMENSION], End[MAX_DIMENSION], Start[MAX_DIMENSION];
  int GridFluxDim[MAX_DIMENSION], GridFluxStartIndex[MAX_DIMENSION];
  float B2;
 
  /* If there are no fields, don't do anything. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifyPhysicalQuantities.\n");
    }
 
    /* If using comoving coordinates, compute a(t) because we'll need it
       to multiply against the CellWidth. */

//    DC revision 16th September 2005 
//    FLOAT a = 1, dadt;
//    if (ComovingCoordinates)
//      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
//        ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
//      }
 
    /* Main loop over all faces. */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* If this dimension is flat, don't do any work. */
 
      if (GridDimension[dim] > 1) {
 
	/* Check that the dims of InitialFluxes & RefinedFluxes are the same */
 
	for (j = 0; j < GridRank; j++)
	  if ((InitialFluxes->LeftFluxStartGlobalIndex[dim][j] !=
	       RefinedFluxes->LeftFluxStartGlobalIndex[dim][j])  ||
	      (InitialFluxes->LeftFluxEndGlobalIndex[dim][j] !=
	       RefinedFluxes->LeftFluxEndGlobalIndex[dim][j])) {
	    ENZO_FAIL("InitialFluxes & RefinedFluxes are different.\n");
	  }
 
	/* Error check Fluxes to make sure they all exist. */
 
	for (field = 0; field < NumberOfBaryonFields; field++)
	  if ((InitialFluxes->LeftFluxes[field][dim] == NULL) ||
	      (RefinedFluxes->LeftFluxes[field][dim] == NULL) ||
	      (InitialFluxes->RightFluxes[field][dim] == NULL) ||
	      (RefinedFluxes->RightFluxes[field][dim] == NULL)) {
	    ENZO_FAIL("Some Flux data is not present.\n");
	  }
 
	/* Compute Start and end indicies of flux region (with respect to
	   the current grid's flux region). */
 
	for (i = 0; i < MAX_DIMENSION; i++) {
	  Start[i] = 0;
	  End[i] = 0;
	}
 
	/* start index = subgrid flux left edge global index -
	                 grid far left edge global index
	     end index = subgrid flux right edge -
	                 grid far left edge global index. */
 
	for (i = 0; i < GridRank; i++) {
	  Start[i] = InitialFluxes->LeftFluxStartGlobalIndex[dim][i] -
	    nlongint((CellLeftEdge[i][0] - DomainLeftEdge[i])/CellWidth[i][0]);
	  End[i] = InitialFluxes->LeftFluxEndGlobalIndex[dim][i] -
	    nlongint((CellLeftEdge[i][0] - DomainLeftEdge[i])/CellWidth[i][0]);
	  if (Start[i] < 0 || End[i] > GridDimension[i]) {
	    fprintf(stderr, "Start/End[%"ISYM"] = %"ISYM"/%"ISYM"\n", dim, Start[i], End[i]);
	    fprintf(stderr, "%"GOUTSYM" %"GOUTSYM" %ld\n",
		    CellLeftEdge[i][0], CellWidth[i][0],
		    InitialFluxes->LeftFluxStartGlobalIndex[dim][i]);
	    ENZO_FAIL("Error in Grid_CorrectForRefinedFluxes!\n");
	  }
	}
 
	/* Correct vector to point at cells just outside the left face.
	   Start[dim] and End[dim] should be the same because the
	   layer to be corrected is but one cell thick. */
 
	Start[dim] = max(Start[dim] - 1, 0);
	End[dim]   = Start[dim];
 
	/* Compute Dimensions of Fluxes */
 
	for (i = 0; i < MAX_DIMENSION; i++)
	  Dim[i] = End[i] - Start[i] + 1;
 
	/* Compute Offset (in baryon field) for right side of face.
	   The +2 is there because we want to correct the cells just the
	   right face.*/
 
	Offset = InitialFluxes->RightFluxStartGlobalIndex[dim][dim] -
	         InitialFluxes->LeftFluxStartGlobalIndex[dim][dim] + 2;
	Offset = min(Offset, GridDimension[dim]-1);  // this isn't needed (?)
 
	/* Check to see if we should correct BoundaryFluxesThisTimeStep
	   instead of the fields themselves. */
 
	CorrectLeftBoundaryFlux = FALSE;
	CorrectRightBoundaryFlux = FALSE;
#ifdef UNUSED
	if (Start[dim] == GridStartIndex[dim]-1)
	  CorrectLeftBoundaryFlux = TRUE;
	if (Start[dim] + Offset == GridEndIndex[dim]+1)
	  CorrectRightBoundaryFlux = TRUE;
#endif /* UNUSED */
 
	/* Set GridFluxStartIndex to the starting position of the flux
	   plane (i.e. exclude grid boundary zones), except for the direction
	   of the flux is set such that GridStartIndex[dim] - Start[dim] = 0 */
	
	for (i = 0; i < MAX_DIMENSION; i++) {
	  GridFluxStartIndex[i] = GridStartIndex[i];
	  GridFluxDim[i] = GridEndIndex[i] - GridStartIndex[i] + 1;
	}
	GridFluxStartIndex[dim] = Start[dim];
	GridFluxDim[dim] = 1;
	
	/* Turn Offset (in dim direction) to Offset (in field array) */
 
	for (i = 0; i < dim; i++)
	  Offset *= GridDimension[i];
 
	/* Multiply faces by density to get conserved quantities
	   (only multiply fields which we are going to correct) */
 
	if (HydroMethod != Zeus_Hydro)
	  for (field = 0; field < NumberOfBaryonFields; field++)
	    if (FieldTypeIsDensity(FieldType[field]) == FALSE &&
		FieldTypeIsRadiation(FieldType[field]) == FALSE &&
		FieldType[field] != Bfield1 &&
		FieldType[field] != Bfield2 && FieldType[field] != Bfield3 &&
		FieldType[field] != PhiField &&
		FieldType[field] != DrivingField1 &&
		FieldType[field] != DrivingField2 &&
		FieldType[field] != DrivingField3 &&
		FieldType[field] != GravPotential)  
	      //		(RadiativeCooling == 0 || (FieldType[field] != TotalEnergy &&
	      //	 			 FieldType[field] != InternalEnergy)))
	      for (k = Start[2]; k <= End[2]; k++)
		for (j = Start[1]; j <= End[1]; j++) {
		  index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
		  for (i = Start[0]; i <= End[0]; i++, index++) {
		    BaryonField[field][index] *= BaryonField[DensNum][index];
		    BaryonField[field][index+Offset] *=
		      BaryonField[DensNum][index+Offset];
		  }
		}
 
	/* Divide species by densities so that at the end we can multiply
           them by the new density (species are not otherwise modified --
	   see the next comment).  This ensures that the species are changed
	   to keep the same fractional density. */
 
	for (field = 0; field < NumberOfBaryonFields; field++)
	  if (FieldType[field] >= ElectronDensity &&
	      FieldType[field] < FieldUndefined &&
	      FieldTypeIsRadiation(FieldType[field]) == FALSE)
	    for (k = Start[2]; k <= End[2]; k++)
	      for (j = Start[1]; j <= End[1]; j++) {
		index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
		for (i = Start[0]; i <= End[0]; i++, index++) {
		  BaryonField[field][index] /= BaryonField[DensNum][index];
		  BaryonField[field][index+Offset] /=
		    BaryonField[DensNum][index+Offset];
		}
	      }
	
	/* Correct face for difference between refined and initial fluxes.
	   (Don't do this for energy if radiative cooling is on because it's
	    no longer conserved.  Similarly, don't do it for the species
	    because they are not individually conserved either -- in fact,
            this could be done for the conserved quantities like charge,
	    total number density summed over ionization, etc.) */
 
	for (field = 0; field < NumberOfBaryonFields; field++)
 	  if (FieldTypeNoInterpolate(FieldType[field]) == FALSE &&
 	      (RadiativeCooling == 0 || (FieldType[field] != TotalEnergy &&
					 FieldType[field] != InternalEnergy))
	      && (FieldType[field] < ElectronDensity)) {
	  for (k = Start[2]; k <= End[2]; k++)
	    for (j = Start[1]; j <= End[1]; j++)
	      for (i = Start[0]; i <= End[0]; i++) {
 
		/* Compute indexes. */
 
		FieldIndex = (k*GridDimension[1] + j)*GridDimension[0] + i;
		FluxIndex  = ((k - Start[2])*Dim[1] + (j - Start[1]))*Dim[0] +
		              (i - Start[0]);
		GridFluxIndex =
                     (i - GridFluxStartIndex[0])
		   + (j - GridFluxStartIndex[1])*GridFluxDim[0]
		   + (k - GridFluxStartIndex[2])*GridFluxDim[1]*GridFluxDim[0];
 
		/* Left side */
 
		if (CorrectLeftBoundaryFlux)
		  /*		  BoundaryFluxesThisTimeStep->LeftFluxes[field][dim][GridFluxIndex] =
		    RefinedFluxes->LeftFluxes[field][dim][FluxIndex]; */
		  ;
		else
		  BaryonField[field][FieldIndex] +=
		     (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		      RefinedFluxes->LeftFluxes[field][dim][FluxIndex] );
 
		if ((FieldTypeIsDensity(FieldType[field]) == TRUE ||
		     FieldType[field] == TotalEnergy ||
		     FieldType[field] == InternalEnergy) &&
		    BaryonField[field][FieldIndex] <= 0) {
		  if (debug)
		    printf("CFRFl warn: %e %e %e %"ISYM" %"ISYM" %"ISYM" %"ISYM" [%"ISYM"]\n",
			   BaryonField[field][FieldIndex],
			   InitialFluxes->LeftFluxes[field][dim][FluxIndex],
			   RefinedFluxes->LeftFluxes[field][dim][FluxIndex],
			   i, j, k, dim, field);
		  
		  /* If new density is < 0 then stop the flux correction. */
 
		  BaryonField[field][FieldIndex] -=
		    (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->LeftFluxes[field][dim][FluxIndex] );

		  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
		    RefinedFluxes->LeftFluxes[ffield][dim][FluxIndex] =
		      InitialFluxes->LeftFluxes[ffield][dim][FluxIndex];
		  
		  ENZO_FAIL("New density or energy is < 0!\n");
		}
 
		/* Right side */
 
		if (CorrectRightBoundaryFlux)
		  /*  BoundaryFluxesThisTimeStep->RightFluxes[field][dim] [GridFluxIndex] =
		  RefinedFluxes->RightFluxes[field][dim][FluxIndex]; */
		  ;
		else
		  BaryonField[field][FieldIndex + Offset] -=
		    (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->RightFluxes[field][dim][FluxIndex] );

 
		if ((FieldTypeIsDensity(FieldType[field]) == TRUE ||
		     FieldType[field] == TotalEnergy ||
		     FieldType[field] == InternalEnergy) &&
		    BaryonField[field][FieldIndex + Offset] <= 0.0) {
		  if (debug)
		    printf("CFRFr warn: %e %e %e %"ISYM" %"ISYM" %"ISYM" %"ISYM" (%"ISYM") [%"ISYM"]\n",
			   BaryonField[field][FieldIndex + Offset],
			   InitialFluxes->RightFluxes[field][dim][FluxIndex],
			   RefinedFluxes->RightFluxes[field][dim][FluxIndex],
			   i, j, k, dim, Offset, field);
 
		  /* If new density is < 0 then stop the flux correction. */
 
		  BaryonField[field][FieldIndex + Offset] +=
		     (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		      RefinedFluxes->RightFluxes[field][dim][FluxIndex] );

		  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
		    RefinedFluxes->RightFluxes[ffield][dim][FluxIndex] =
		      InitialFluxes->RightFluxes[ffield][dim][FluxIndex];

		  ENZO_FAIL("New density or energy is < 0!\n");
		}
 
	      }
	  }
 
	/* Return faces to original quantity. */
 
	if (HydroMethod != Zeus_Hydro)
	  for (field = 0; field < NumberOfBaryonFields; field++)
	    if (FieldTypeIsDensity(FieldType[field]) == FALSE &&
		FieldTypeIsRadiation(FieldType[field]) == FALSE &&
		FieldType[field] != Bfield1 &&
		FieldType[field] != Bfield2 && FieldType[field] != Bfield3 &&
		FieldType[field] != PhiField &&
		FieldType[field] != DrivingField1 &&
		FieldType[field] != DrivingField2 &&
		FieldType[field] != DrivingField3 &&
		FieldType[field] != GravPotential)
	      //		(RadiativeCooling == 0 || (FieldType[field] != TotalEnergy &&
	      //	 			 FieldType[field] != InternalEnergy)))
	      for (k = Start[2]; k <= End[2]; k++)
		for (j = Start[1]; j <= End[1]; j++) {
		  index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
		  for (i = Start[0]; i <= End[0]; i++, index++) {
		    BaryonField[field][index] /= BaryonField[DensNum][index];
		    BaryonField[field][index+Offset] /=
		      BaryonField[DensNum][index+Offset];
		  }
		}
	 
	/* If appropriate, restore consistency between total and internal
	   energy in corrected faces. */
 
	if (DualEnergyFormalism == TRUE && RadiativeCooling == FALSE) {
	  for (k = Start[2]; k <= End[2]; k++)
	    for (j = Start[1]; j <= End[1]; j++) {
	      i1 = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
	      i2 = i1 + Offset;
	      for (i = Start[0]; i <= End[0]; i++, i1++, i2++) {
		BaryonField[GENum][i1] = max(BaryonField[GENum][i1],
					     tiny_number);
		BaryonField[GENum][i2] = max(BaryonField[GENum][i2],
					     tiny_number);
		BaryonField[TENum][i1] = BaryonField[GENum][i1] +
		  0.5 * BaryonField[Vel1Num][i1] * BaryonField[Vel1Num][i1];
		BaryonField[TENum][i2] = BaryonField[GENum][i2] +
		  0.5 * BaryonField[Vel1Num][i2] * BaryonField[Vel1Num][i2];
		if (GridRank > 1) {
		  BaryonField[TENum][i1] +=
		    0.5 * BaryonField[Vel2Num][i1] * BaryonField[Vel2Num][i1];
		  BaryonField[TENum][i2] +=
		    0.5 * BaryonField[Vel2Num][i2] * BaryonField[Vel2Num][i2];
		}
		if (GridRank > 2) {
		  BaryonField[TENum][i1] +=
		    0.5 * BaryonField[Vel3Num][i1] * BaryonField[Vel3Num][i1];
		  BaryonField[TENum][i2] +=
		    0.5 * BaryonField[Vel3Num][i2] * BaryonField[Vel3Num][i2];
		}
		if (HydroMethod == MHD_RK) {
		  B2 = POW(BaryonField[B1Num][i1],2) + POW(BaryonField[B2Num][i1],2) +
		    POW(BaryonField[B3Num][i1],2);
		  BaryonField[TENum][i1] += 0.5 * B2 / BaryonField[DensNum][i1];
		  B2 = POW(BaryonField[B1Num][i2],2) + POW(BaryonField[B2Num][i2],2) +
		    pow(BaryonField[B3Num][i2],2);
		  BaryonField[TENum][i2] += 0.5 * B2 / BaryonField[DensNum][i2];
		}

 
	      }
		
	    } // end: loop over faces
 
	} // end: if (DualEnergyFormalism)
 
	/* Multiply species by density to return from fractional to real
	   density. (see comments above regarding species). */
 
	for (field = 0; field < NumberOfBaryonFields; field++)
	  if (FieldType[field] >= ElectronDensity &&

	      FieldType[field] < FieldUndefined &&
	      FieldTypeIsRadiation(FieldType[field]) == FALSE)
	    for (k = Start[2]; k <= End[2]; k++)
	      for (j = Start[1]; j <= End[1]; j++) {
		index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
		for (i = Start[0]; i <= End[0]; i++, index++) {
		  BaryonField[field][index] *= BaryonField[DensNum][index];
		  BaryonField[field][index+Offset] *=
		    BaryonField[DensNum][index+Offset];
		}
	      }
 
      }  // end: if GridDimension[dim] > 1
 
      /* delete Refined fluxes as they're not needed anymore. */
 
      for (field = 0; field < NumberOfBaryonFields; field++) {
	delete RefinedFluxes->LeftFluxes[field][dim];
	delete RefinedFluxes->RightFluxes[field][dim];
	RefinedFluxes->LeftFluxes[field][dim] = NULL;
	RefinedFluxes->RightFluxes[field][dim] = NULL;
      }
 
    } // next dimension
 
  }
 
  return SUCCESS;
 
}
#endif
