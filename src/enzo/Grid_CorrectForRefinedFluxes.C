/***********************************************************************
/
/  GRID CLASS (CORRECT SOLUTION GIVEN ORIGINAL AND REFINED FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       January, 2003
/	       Include extra fields beyond Metallicity!
/  modified2:  David Collins & Rick Wagner
/  date:       May, 2005
/	       Include flux correction for outside grids.
/              Re-instated CorrectBoundaryFlux code.
/  modified2: David Collins, 2005
/              Updated algebra so Cosmological Expansion is also
/              conservative.  This fix also came with fixes to euler.src and
/              Grid_GetProjectedBoundaryFluxes.C, so make sure you get those.
/
/  PURPOSE:    Ensures conservation of stuff.
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
//   Also note that subgrids sharing a face with This Grid, but not proper subgrids,
//   also need to be taken into account.
 
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
int MakeFieldConservative(field_type field); 
 
int grid::CorrectForRefinedFluxes(fluxes *InitialFluxes,
				  fluxes *RefinedFluxes,
				  fluxes *BoundaryFluxesThisTimeStep,
				  int SUBlingGrid,
				  TopGridData *MetaData)
{
 
  // Return if this doesn't concern us.
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro)
    return SUCCESS;
 
  // declarations
 
  int i1, i2, i, j, k, dim, field, ffield, index;
  int FieldIndex, FluxIndex, GridFluxIndex, Offset, RefinedFluxIndex;
  int End[MAX_DIMENSION], Start[MAX_DIMENSION];
 
  //Dimensions of Initial and Refined fluxes.  ("Dim" should be "InitialDim")
  int Dim[MAX_DIMENSION],RefinedDim[MAX_DIMENSION] = {0,0,0};
 
  //Dimension and Start Index for the ParentGrid (this grid) boundary flux.
  int GridFluxDim[MAX_DIMENSION], GridFluxStartIndex[MAX_DIMENSION];
 
  // Used for calculating position in the RefinedFlux and InitialFlux data structure.
  int RefinedOffset[MAX_DIMENSION] ={0,0,0}, InitialOffset[MAX_DIMENSION] = {0,0,0};
 
  //Internal flags: Correct the BaryonField
  int CorrectLeftBaryonField, CorrectRightBaryonField;
  //For correction of the Parent Grid (this grid) boundary flux.
  int CorrectLeftBoundaryFlux, CorrectRightBoundaryFlux;
  float CorrectionAmountLeft, CorrectionAmountRight; 
 
  long_int GlobalDim;
 
  /* If there are no fields, don't do anything. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifyPhysicalQuantities.");
    }

    //dcc kludge:  Just remove a(t)? 09/06/05 
    /* If using comoving coordinates, compute a(t) because we'll need it
       to multiply against the CellWidth. */

    //  DC revision September 16th 2005 
    //    FLOAT a = 1, dadt;
    //    if (ComovingCoordinates)
    //      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    //        ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    //      }
 
 
    /* Main loop over all faces. */
 
    for (dim = 0; dim < GridRank; dim++) {
      if (GridDimension[dim] > 1) {
	
	/* Check that the dims of InitialFluxes & RefinedFluxes are the same */
 
	/* don't do this for SUBlings */
	if( SUBlingGrid == FALSE ){
	  for (j = 0; j < GridRank; j++)
	    if ((InitialFluxes->LeftFluxStartGlobalIndex[dim][j] !=
		 RefinedFluxes->LeftFluxStartGlobalIndex[dim][j])  ||
		(InitialFluxes->LeftFluxEndGlobalIndex[dim][j] !=
		 RefinedFluxes->LeftFluxEndGlobalIndex[dim][j])) {
//	      printf("dim=%d / j=%d //// %d == %d :: %d == %d\n", 
//		     dim, j, 
//		     InitialFluxes->LeftFluxStartGlobalIndex[dim][j],
//		     RefinedFluxes->LeftFluxStartGlobalIndex[dim][j],
//		     InitialFluxes->LeftFluxEndGlobalIndex[dim][j],
//		     RefinedFluxes->LeftFluxEndGlobalIndex[dim][j]);
	      ENZO_FAIL("InitialFluxes & RefinedFluxes are different.\n");
	    }
	  /* Error check Fluxes to make sure they all exist. */
	  for (field = 0; field < NumberOfBaryonFields; field++)
	    if ((InitialFluxes->LeftFluxes[field][dim] == NULL) ||
		(RefinedFluxes->LeftFluxes[field][dim] == NULL) ||
		(InitialFluxes->RightFluxes[field][dim] == NULL) ||
		(RefinedFluxes->RightFluxes[field][dim] == NULL)) {
	      fprintf(stderr,"Some Flux data is not present.\n");
	    return FAIL;
	  }
	}
 
	//by default, we want to correct the flux.
	CorrectLeftBaryonField = CorrectRightBaryonField = TRUE;
	if( SUBlingGrid == TRUE ){
	
	  /* calculate Global dimensions on this level */
	  GlobalDim = 	nlongint(( DomainRightEdge[dim] - DomainLeftEdge[dim])
				 / CellWidth[dim][0]);
	
	  /* get the dims of the refined fluxes to calculate
	     array indices */
	
	  for (i = 0; i < MAX_DIMENSION; i++){
	    RefinedDim[i] = RefinedFluxes->LeftFluxEndGlobalIndex[dim][i] -
	      RefinedFluxes->LeftFluxStartGlobalIndex[dim][i] + 1;
	  }
	
	  /* check if SUBling left or right edge lies on
	     Domain boundary, and if so, modulo the indices,
	     otherwise, bump the indices to match the initial flux's.
	  */
	
	  if( RefinedFluxes->LeftFluxStartGlobalIndex[dim][dim] == 0 &&
	      MetaData->LeftFaceBoundaryCondition[dim] == periodic ){
	    RefinedFluxes->LeftFluxStartGlobalIndex[dim][dim] = GlobalDim - 1;
	    RefinedFluxes->LeftFluxEndGlobalIndex[dim][dim] = GlobalDim - 1;
	  }else{
	    RefinedFluxes->LeftFluxStartGlobalIndex[dim][dim]--;
	    RefinedFluxes->LeftFluxEndGlobalIndex[dim][dim]--;
	  }
	
	
	  if( RefinedFluxes->RightFluxStartGlobalIndex[dim][dim] == GlobalDim - 1 &&
	      MetaData->RightFaceBoundaryCondition[dim] == periodic){
	    RefinedFluxes->RightFluxStartGlobalIndex[dim][dim] = 0;
	    RefinedFluxes->RightFluxEndGlobalIndex[dim][dim] = 0;
	  }else{
	    RefinedFluxes->RightFluxStartGlobalIndex[dim][dim]++;
	    RefinedFluxes->RightFluxEndGlobalIndex[dim][dim]++;
	  }
	
	  /* check to see if we're doing this dimension at all.
	     only the dimension of contact needs to be checked,
	     since SUBling grids can only have contact along a
	     single axis. any corrections to this statement
	     earns a beer */
	
	
	  //note these are integers, so comparing them directly is ok.
	  //Also note that here we don't do the complete check for SUBling-ness.
	  //More logic is necessary for any two arbitrary grids, but it's been done
	  //already in populating the SUBling list.  Here, all we need
	  //to do is determine which face needs the correction, as we already know one exists.
	
	  if( InitialFluxes->RightFluxStartGlobalIndex[dim][dim] ==
	      RefinedFluxes->LeftFluxStartGlobalIndex[dim][dim] ){
	    CorrectRightBaryonField = TRUE;
	  }else{
	    CorrectRightBaryonField = FALSE;
	  }
	
	  if( InitialFluxes->LeftFluxStartGlobalIndex[dim][dim]==
	      RefinedFluxes->RightFluxStartGlobalIndex[dim][dim] ){
	    CorrectLeftBaryonField = TRUE;
	  }else{
	    CorrectLeftBaryonField = FALSE;
	  }
	
	  for (i = 0; i < MAX_DIMENSION; i++) {
	    /* calculate the offset, so the index of the refined fluxes can
	       be determined from the grid's index */
	    RefinedOffset[i] = max(InitialFluxes->LeftFluxStartGlobalIndex[dim][i]-
				   RefinedFluxes->LeftFluxStartGlobalIndex[dim][i],0);
	  }
 
	  RefinedFluxes->LeftFluxStartGlobalIndex[dim][dim]=0;
	  RefinedOffset[dim]=0;
 
	}//Subling == TRUE	
 
	if( CorrectLeftBaryonField || CorrectRightBaryonField){
	
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
	
	  /* modified to account for different dimensions to the
	     initial and refined fluxes */
	
	  for (i = 0; i < GridRank; i++) {
	    Start[i] = max(InitialFluxes->LeftFluxStartGlobalIndex[dim][i],
			   RefinedFluxes->LeftFluxStartGlobalIndex[dim][i]) -
	      nlongint((CellLeftEdge[i][0] - DomainLeftEdge[i])/
		       CellWidth[i][0]);
	    End[i] = min(InitialFluxes->LeftFluxEndGlobalIndex[dim][i],
			 RefinedFluxes->LeftFluxEndGlobalIndex[dim][i]) -
	      nlongint((CellLeftEdge[i][0] - DomainLeftEdge[i])/
		       CellWidth[i][0]);
	
	
	    if (Start[i] < 0 || End[i] > GridDimension[i]) {
	      fprintf(stderr, "Start/End[%"ISYM"] = %"ISYM"/%"ISYM"\n",
		      dim, Start[i], End[i]);
	      fprintf(stderr, "%"GOUTSYM" %"GOUTSYM" %lld\n",
		      CellLeftEdge[i][0], CellWidth[i][0],
		      InitialFluxes->LeftFluxStartGlobalIndex[dim][i]);
	      ENZO_FAIL("Error in FluxFix_Grid_CorrectForRefinedFluxes!\n");
	    }
	  }
	
	  /* Correct vector to point at cells just outside the left face.
	     Start[dim] and End[dim] should be the same because the
	     layer to be corrected is but one cell thick. */
	
	  Start[dim] = max(Start[dim] - 1, 0);
	  End[dim]   = Start[dim];
 
	  /* Compute Dimensions of InitialFluxes */
	
	  for (i = 0; i < MAX_DIMENSION; i++)
	    Dim[i] = End[i] - Start[i] + 1;
	
	  /* Compute Offset (in baryon field) for right side of face.
	     The +2 is there because we want to correct the cells just the
	     right face.*/
	
	  Offset = InitialFluxes->RightFluxStartGlobalIndex[dim][dim] -
	    InitialFluxes->LeftFluxStartGlobalIndex[dim][dim] + 2;
	  Offset = min(Offset, GridDimension[dim]-1);  // this isn't needed (?)
 
 
	  //For SUBling grids, alter Offset, Start, and End to reflect that we're
	  //adjusting the INNER edge of the grid if the SUBgrid is outside of it.
 
	  if( SUBlingGrid ){
	    Offset -= 2;
	    Start[dim]++;
	    End[dim]++;
	
	    //Also correct Dim, the size of the Initial Flux: it comes from This Grid,
	    //not the Subgrid.
	
	    for(i=0;i<GridRank;i++)
	      if(i != dim){
		Dim[i] = GridEndIndex[i]-GridStartIndex[i]+1;
		InitialOffset[i] = max( RefinedFluxes->LeftFluxStartGlobalIndex[dim][i]-
					InitialFluxes->LeftFluxStartGlobalIndex[dim][i],
					0);
	      }else{
		Dim[i] = 1;
		InitialOffset[i] = 0;
	      }
	  }
	
	
	  /* Check to see if we should correct BoundaryFluxesThisTimeStep
	     instead of the fields themselves. */
	
	  CorrectLeftBoundaryFlux = FALSE;
	  CorrectRightBoundaryFlux = FALSE;
	
	  if (Start[dim] == GridStartIndex[dim]-1){
	    CorrectLeftBoundaryFlux = TRUE;
	    CorrectLeftBaryonField  = FALSE;
	  }
	  if (Start[dim] + Offset == GridEndIndex[dim]+1){
	    CorrectRightBoundaryFlux = TRUE;
	    CorrectRightBaryonField  = FALSE;
	  }
 
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
	
	  if (HydroMethod != Zeus_Hydro) {

	    for (field = 0; field < NumberOfBaryonFields; field++) 
	      if ( MakeFieldConservative(FieldType[field]) ) {
		for (k = Start[2]; k <= End[2]; k++) {
		  for (j = Start[1]; j <= End[1]; j++) {
		    index = (k*GridDimension[1] + j)*GridDimension[0] + Start[0];
		    for (i = Start[0]; i <= End[0]; i++, index++) {
		      BaryonField[field][index] *= BaryonField[DensNum][index];
		      BaryonField[field][index+Offset] *= 
			BaryonField[DensNum][index+Offset];
		    }
		  }
		}
	      }
	  }

	  
	
	  /* Divide species by densities so that at the end we can multiply
	     them by the new density (species are not otherwise modified --
	     see the next comment).  This ensures that the species are changed
	     to keep the same fractional density. */
	
	  for (field = 0; field < NumberOfBaryonFields; field++)
	    if (FieldType[field] >= ElectronDensity &&
		FieldType[field] <= Metallicity &&
		FieldTypeNoInterpolate(FieldType[field]) == FALSE &&
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
	
	  if (Coordinate == Cartesian) {
	  for (field = 0; field < NumberOfBaryonFields; field++){
	    if ((FieldTypeNoInterpolate(FieldType[field]) == FALSE) &&
		(RadiativeCooling == 0 || (FieldType[field] != TotalEnergy &&
					   FieldType[field] != InternalEnergy)) &&
		(FieldType[field] < ElectronDensity) && 
		FieldType[field] != DrivingField1 &&
		FieldType[field] != DrivingField2 &&
		FieldType[field] != DrivingField3 &&
		FieldType[field] != GravPotential &&
		FieldType[field] != DebugField) {
	      for (k = Start[2]; k <= End[2]; k++){
		for (j = Start[1]; j <= End[1]; j++){
		  for (i = Start[0]; i <= End[0]; i++) {
		
		    /* Compute indexes. */
		
		    FieldIndex = (k*GridDimension[1] + j)*GridDimension[0] + i;
		    FluxIndex  = ((k - Start[2]+InitialOffset[2])*Dim[1] +
				  (j - Start[1]+InitialOffset[1]))*Dim[0] +
		      (i - Start[0]+InitialOffset[0]);
		
		    if( SUBlingGrid ){
		      RefinedFluxIndex = ((k - Start[2] + RefinedOffset[2])*RefinedDim[1] +
					  (j - Start[1] + RefinedOffset[1]))*RefinedDim[0] +
			(i - Start[0] + RefinedOffset[0]);
		
		    }else{
		      RefinedFluxIndex = FluxIndex;
		    }
		
		    GridFluxIndex =
		      (i - GridFluxStartIndex[0])
		      + (j - GridFluxStartIndex[1])*GridFluxDim[0]
		      + (k - GridFluxStartIndex[2])*GridFluxDim[1]*GridFluxDim[0];
		
		
 		    if (CorrectLeftBoundaryFlux && CorrectParentBoundaryFlux )
		      BoundaryFluxesThisTimeStep->LeftFluxes[field][dim][GridFluxIndex] =
			RefinedFluxes->LeftFluxes[field][dim][FluxIndex];
		
		    if(CorrectLeftBaryonField){
		
		      if( SUBlingGrid == FALSE ){
			CorrectionAmountLeft = 
			  InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
			  RefinedFluxes->LeftFluxes[field][dim][FluxIndex];
			BaryonField[field][FieldIndex] += CorrectionAmountLeft;
			
		      }else{ /* if( SUBlingGrid == False) */
			
			CorrectionAmountLeft = 
			  -(InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
			    RefinedFluxes->RightFluxes[field][dim][RefinedFluxIndex]);
			BaryonField[field][FieldIndex] += CorrectionAmountLeft;
			
		      }
		    }
		
 		    if (CorrectRightBoundaryFlux && CorrectParentBoundaryFlux )
		      BoundaryFluxesThisTimeStep->RightFluxes[field][dim] [GridFluxIndex] =
			RefinedFluxes->RightFluxes[field][dim][FluxIndex];
		
		    /* update only if necessary */
		    if(CorrectRightBaryonField){
		
		      if( SUBlingGrid == FALSE ){

			CorrectionAmountRight = 
			  -(InitialFluxes->RightFluxes[field][dim][FluxIndex] -
			    RefinedFluxes->RightFluxes[field][dim][FluxIndex]);
			  
			BaryonField[field][FieldIndex + Offset] += CorrectionAmountRight;
			
		      }else{ /* if( SUBlingGrid == FALSE ){ */
			  CorrectionAmountRight = 
			    InitialFluxes->RightFluxes[field][dim][FluxIndex] -
			    RefinedFluxes->LeftFluxes[field][dim][RefinedFluxIndex];
			  BaryonField[field][FieldIndex + Offset] += CorrectionAmountRight;
			
		      } // else{ /* if( SUBlingGrid == FALSE ){ */
		    } // if(CorrectRightBaryonField)
		
		    if ((FieldTypeIsDensity(FieldType[field]) == TRUE ||
			 FieldType[field] == TotalEnergy ||
			 FieldType[field] == InternalEnergy)) {

		      /* If new density & energy is < 0 then undo the
			 flux correction. */

		      if (CorrectLeftBaryonField &&
			  BaryonField[field][FieldIndex] <= 0) {
			BaryonField[field][FieldIndex] -= CorrectionAmountLeft;

			if (SUBlingGrid == FALSE) {
			  if (debug)
			    printf("P(%d) -- CFRFl warn: %e %e %e %e %"ISYM
				   " %"ISYM" %"ISYM" %"ISYM" [%"ISYM"]\n",
				   MyProcessorNumber, BaryonField[field][FieldIndex],
				   InitialFluxes->LeftFluxes[field][dim][FluxIndex],
				   RefinedFluxes->LeftFluxes[field][dim][FluxIndex],
				   CorrectionAmountLeft,
				   i, j, k, dim, field);
			  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
			    RefinedFluxes->LeftFluxes[ffield][dim][FluxIndex] =
			      InitialFluxes->LeftFluxes[ffield][dim][FluxIndex];
			} else {
			  if (debug)
			    printf("P(%d) -- CFRFlS warn: %e %e %e %e %"ISYM
				   " %"ISYM" %"ISYM" %"ISYM" [%"ISYM"]\n",
				   MyProcessorNumber, BaryonField[field][FieldIndex],
				   InitialFluxes->LeftFluxes[field][dim][FluxIndex],
				   RefinedFluxes->RightFluxes[field][dim][RefinedFluxIndex],
				   CorrectionAmountLeft,
				   i, j, k, dim, field);
			  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
			    RefinedFluxes->RightFluxes[ffield][dim][RefinedFluxIndex] =
			      InitialFluxes->LeftFluxes[ffield][dim][FluxIndex];
			} // ENDELSE (SUBlingGrid == FALSE)
		      } // ENDIF CorrectLeftBaryonField

		      if (CorrectRightBaryonField &&
			  BaryonField[field][FieldIndex+Offset] <= 0) {
			BaryonField[field][FieldIndex+Offset] -= CorrectionAmountRight;

			if (SUBlingGrid == FALSE) {
			  if (debug)
			    printf("P(%d) -- CFRFr warn: %e %e %e %e %"ISYM
				   " %"ISYM" %"ISYM" %"ISYM" [%"ISYM"]\n",
				   MyProcessorNumber, BaryonField[field][FieldIndex],
				   InitialFluxes->RightFluxes[field][dim][FluxIndex],
				   RefinedFluxes->RightFluxes[field][dim][FluxIndex],
				   CorrectionAmountRight,
				   i, j, k, dim, field);
			  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
			    RefinedFluxes->RightFluxes[ffield][dim][FluxIndex] =
			      InitialFluxes->RightFluxes[ffield][dim][FluxIndex];
			} else {
			  if (debug)
			    printf("P(%d) -- CFRFrS warn: %e %e %e %e %"ISYM
				   " %"ISYM" %"ISYM" %"ISYM" [%"ISYM"]\n",
				   MyProcessorNumber, BaryonField[field][FieldIndex],
				   InitialFluxes->LeftFluxes[field][dim][FluxIndex],
				   RefinedFluxes->RightFluxes[field][dim][RefinedFluxIndex],
				   CorrectionAmountRight,
				   i, j, k, dim, field);
			  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
			    RefinedFluxes->LeftFluxes[ffield][dim][RefinedFluxIndex] =
			      InitialFluxes->RightFluxes[ffield][dim][FluxIndex];
			}
		      } // ENDIF CorrectRightBaryonField
		    }
		  }// for (i = Start[0]; i <= End[0]; i++) {
		} // for (j = Start[1]; j <= End[1]; j++){
	      } // for (k = Start[2]; k <= End[2]; k++){
	    }	// if ((RadiativeCooling == 0 || (FieldType[field] != TotalEnergy && etc
	  } // for (field = 0; field < NumberOfBaryonFields; field++){
	  }


	if (Coordinate == Cylindrical) {
	FLOAT xr, xl, xc, geofacr, geofacl;
	for (field = 0; field < NumberOfBaryonFields; field++) {
	  /*if ((RadiativeCooling == 0 || (FieldType[field] != TotalEnergy && 
	    FieldType[field] != InternalEnergy))
	    && (FieldType[field] < ElectronDensity) ) {*/
	  if (FieldType[field] < ElectronDensity &&
	      FieldType[field] != DrivingField1 &&
	      FieldType[field] != DrivingField2 &&
	      FieldType[field] != DrivingField3 &&
	      FieldType[field] != GravPotential &&
	      FieldType[field] != DebugField
	      ) {
	  for (k = Start[2]; k <= End[2]; k++) {
	    for (j = Start[1]; j <= End[1]; j++) {
	      for (i = Start[0]; i <= End[0]; i++) {
		/* Compute indexes. */
		FieldIndex = (k*GridDimension[1] + j)*GridDimension[0] + i;
		FluxIndex  = ((k - Start[2])*Dim[1] + (j - Start[1]))*Dim[0] +
		              (i - Start[0]);
		GridFluxIndex = 
                     (i - GridFluxStartIndex[0]) 
		   + (j - GridFluxStartIndex[1])*GridFluxDim[0]
		   + (k - GridFluxStartIndex[2])*GridFluxDim[1]*GridFluxDim[0];
		if (dim == 0) {
		  /* Left side */
		  xr = CellLeftEdge[0][i] + CellWidth[0][i];
                  xc = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
                  geofacr = xr/xc;
		  BaryonField[field][FieldIndex] +=
		    (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->LeftFluxes[field][dim][FluxIndex] )*geofacr;

		  /* Right side */
		  xl = CellLeftEdge[0][i+Offset];
                  xc = xl + 0.5*CellWidth[0][i+Offset];
                  geofacl = xl/xc;
		  BaryonField[field][FieldIndex + Offset] -=
		    (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->RightFluxes[field][dim][FluxIndex] )*geofacl;
		}		
		if (dim == 1) {
		  /* Left side */
		  BaryonField[field][FieldIndex] +=
		    (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->LeftFluxes[field][dim][FluxIndex] );
		  /* Right side */
		  BaryonField[field][FieldIndex + Offset] -=
		    (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->RightFluxes[field][dim][FluxIndex] );
		}		
		if (dim == 2) {

		  /* Left side */
		  xc = CellLeftEdge[0][i] + 0.5*CellWidth[0][i]; 
		  BaryonField[field][FieldIndex] +=
		    (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->LeftFluxes[field][dim][FluxIndex] )*xc;
		  /* Right side */
		  xl = CellLeftEdge[0][i+Offset];
                  xc = xl + 0.5*CellWidth[0][i+Offset];
		  BaryonField[field][FieldIndex + Offset] -=
		    (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		     RefinedFluxes->RightFluxes[field][dim][FluxIndex] )*xc;

		}		
		/* Check for posivtivity and undo flux correction if negative */
		if ((FieldType[field] == Density || 
		     FieldType[field] == TotalEnergy ||
		     FieldType[field] == InternalEnergy) &&
		    BaryonField[field][FieldIndex] <= 0) {
		  /*if (debug) {
		    printf("CFRFl warn: %e %e %e %d %d %d %d [%d]\n",
			   BaryonField[field][FieldIndex],
			   InitialFluxes->LeftFluxes[field][dim][FluxIndex],
			   RefinedFluxes->LeftFluxes[field][dim][FluxIndex],
			   i, j, k, dim, field);
			   }*/
		  /* If new density & energy is < 0 then undo the flux correction. */
		  BaryonField[field][FieldIndex] -=
		     (InitialFluxes->LeftFluxes[field][dim][FluxIndex] -
		      RefinedFluxes->LeftFluxes[field][dim][FluxIndex] );
		  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
		    RefinedFluxes->LeftFluxes[ffield][dim][FluxIndex] =
		      InitialFluxes->LeftFluxes[ffield][dim][FluxIndex];
		}
		if ((FieldType[field] == Density || 
		     FieldType[field] == TotalEnergy ||
		     FieldType[field] == InternalEnergy) &&
		    BaryonField[field][FieldIndex + Offset] <= 0.0) {
		  /*if (debug) {
		    printf("CFRFr warn: %e %e %e %d %d %d %d (%d) [%d]\n",
			   BaryonField[field][FieldIndex + Offset],
			   InitialFluxes->RightFluxes[field][dim][FluxIndex],
			   RefinedFluxes->RightFluxes[field][dim][FluxIndex],
			   i, j, k, dim, Offset, field);
			   }*/
		  /* If new density & energy is < 0 then undo the flux correction. */
		  BaryonField[field][FieldIndex + Offset] +=
		     (InitialFluxes->RightFluxes[field][dim][FluxIndex] -
		      RefinedFluxes->RightFluxes[field][dim][FluxIndex] );
		  for (ffield = 0; ffield < NumberOfBaryonFields; ffield++)
		    RefinedFluxes->RightFluxes[ffield][dim][FluxIndex] =
		      InitialFluxes->RightFluxes[ffield][dim][FluxIndex];
		}
	      }
	    }
	  }
	  
	  }
	}
	}
	

	    /* Return faces to original quantity. */
	
	  if (HydroMethod != Zeus_Hydro)
	    for (field = 0; field < NumberOfBaryonFields; field++)
	      if ( MakeFieldConservative(FieldType[field]))
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
	
	  if (DualEnergyFormalism == TRUE && RadiativeCooling == FALSE){
	    float B2 = 0.0;
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
		    B2 = POW(BaryonField[B1Num][i1],2) + 
		      POW(BaryonField[B2Num][i1],2) +
		      POW(BaryonField[B3Num][i1],2);
		    BaryonField[TENum][i1] += 0.5 * B2 / BaryonField[DensNum][i1];
		    B2 = POW(BaryonField[B1Num][i2],2) + 
		      POW(BaryonField[B2Num][i2],2) +
		      POW(BaryonField[B3Num][i2],2);
		    BaryonField[TENum][i2] += 0.5 * B2 / BaryonField[DensNum][i2];
		  }

		  if (HydroMethod == MHD_Li){
		    B2 = POW(CenteredB[0][i1],2) + 
		      POW(CenteredB[1][i1],2) +
		      POW(CenteredB[2][i1],2);
		    BaryonField[TENum][i1] += 0.5 * B2 / BaryonField[DensNum][i1];
		    B2 = POW(CenteredB[0][i2],2) + 
		      POW(CenteredB[1][i2],2) +
		      POW(CenteredB[2][i2],2);
		    BaryonField[TENum][i2] += 0.5 * B2 / BaryonField[DensNum][i2];
		  }
		
		}		
	      } // end: loop over faces
	  } // end: if (DualEnergyFormalism)
	
	    /* Multiply species by density to return from fractional to real
	       density. (see comments above regarding species). */
	
	  for (field = 0; field < NumberOfBaryonFields; field++)
	    if (FieldType[field] >= ElectronDensity &&
		FieldType[field] <= Metallicity &&
		FieldTypeNoInterpolate(FieldType[field]) == FALSE &&
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
	
	} // if( CorrectLeftBaryonField || CorrectRightBaryonField){

      } // end: if GridDimension[dim] > 1
      /* delete Refined fluxes as they're not needed anymore. */
 
      for (field = 0; field < NumberOfBaryonFields; field++) {
	delete [] RefinedFluxes->LeftFluxes[field][dim];
	delete [] RefinedFluxes->RightFluxes[field][dim];
	RefinedFluxes->LeftFluxes[field][dim] = NULL;
	RefinedFluxes->RightFluxes[field][dim] = NULL;
      }
 
    } // next dimension
  } // Number of baryons fields > 0
 
  return SUCCESS;
 
}
