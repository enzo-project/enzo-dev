#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber
/***********************************************************************
/
/  GRID CLASS (TRANSPORT PHOTON PACKAGES)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This is the heart of the radiative transfer algorithm.
/    On each Grid we initialize photo and heating rates and then call
/    WalkPhotonPackage so all photon packages are transported along their
/    own directions and the photo-ionization and heating rates on 
/    on the grid are updated on the fly. 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "phys_constants.h"

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
PhotonPackageEntry *PopPhoton(PhotonPackageEntry * &Node);
PhotonPackageEntry *DeletePhotonPackage(PhotonPackageEntry *PP);
int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::TransportPhotonPackages(int level, int finest_level, 
				  ListOfPhotonsToMove **PhotonsToMove, 
				  int GridNum, grid **Grids0, int nGrids0, 
				  grid *ParentGrid, grid *CurrentGrid)
{

  int i,j,k, dim, index, count;
  grid *MoveToGrid;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || MultiSpecies < 1) 
    return SUCCESS;

  if (RadiativeTransfer < 1) 
    return SUCCESS;

  if (RadiativeTransfer > 0 && GridRank < 3) {
    ENZO_FAIL("Transfer in less than 3D is not implemented!\n");
  }

  if (PhotonPackages->NextPackage == NULL)
    return SUCCESS;

  /* Get units. */
  double MassUnits, RT_Units;
  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  MassUnits = (double) DensityUnits * POW(LengthUnits, 3.0);
  RT_Units = (double) TimeUnits * POW(LengthUnits, -3.0);

  /* speed of light in code units. note this one is independent of
     a(t), and Modify the photon propagation speed by this
     parameter */

  const float c_cgs = 2.99792e10;
  float LightSpeed;
  LightSpeed = RadiativeTransferPropagationSpeedFraction * 
    (c_cgs/VelocityUnits);

  float DomainWidth[MAX_DIMENSION];
  //double FinestCellVolume = pow(RefineBy, -3*(finest_level-level));
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    //FinestCellVolume *= CellWidth[dim][0];
  }

  // if (DEBUG) fprintf(stdout,"TransportPhotonPackage: initialize fields.\n");
  // if (DEBUG) fprintf(stdout,"TransportPhotonPackage: %"ISYM" %"ISYM" .\n",
  // 		     GridStartIndex[0], GridEndIndex[0]);

  PhotonPackageEntry *PP, *FPP, *SavedPP, *PausedPP;
  PP = PhotonPackages;

  if (DEBUG) {
    count = 0;
    while ((PP->NextPackage) != NULL) { 
      count++;
      PP=PP->NextPackage;
    }
    fprintf(stdout, "TransportPhotonPackage: done initializing.\n");
    fprintf(stdout, "[%d] counted %"ISYM" packages\n", this->ID, count);
  }

  /* If requested, make vertex centered field (only when it doesn't
     exist ... see inside routine). */

  if (RadiativeTransferInterpolateField)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (FieldsToInterpolate[i] == TRUE)
	if (this->ComputeVertexCenteredField(i) == FAIL) {
	  ENZO_VFAIL("Error in grid->ComputeVertexCenteredField "
		  "(field %"ISYM").\n", i)
	}

  /* Calculate minimum photon flux before a ray is deleted */
  
  const double alphaB = 2.6e-13;
  float MinimumPhotonFlux, RecombinationTime;
  int gmethod = INT_UNDEFINED;
  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
    if (CellFlaggingMethod[i] == 2) gmethod = i;
  if (gmethod == INT_UNDEFINED)
    MinimumPhotonFlux = POW(TopGridDx[0], GridRank) * POW(RefineBy, level); // estimate
  else
    MinimumPhotonFlux = MinimumMassForRefinement[gmethod] * 
      POW(RefineBy, level*MinimumMassForRefinementLevelExponent[gmethod]);
  if (ComovingCoordinates)
    MinimumPhotonFlux *= (float) ((RT_Units / TimeUnits) * (MassUnits / mh) * dtPhoton / 
				  (PhotonTime * RadiativeTransferHubbleTimeFraction));
  else {
    RecombinationTime = 1.0 / (alphaB * DensityUnits / mh) / TimeUnits;
    MinimumPhotonFlux *= (float) ((RT_Units / TimeUnits) * (MassUnits / mh) * dtPhoton / 
				  (10*RecombinationTime));
  }
  
  count = 0;
  PP = PhotonPackages->NextPackage;
  FPP = this->FinishedPhotonPackages;
  PausedPP = this->PausedPhotonPackages;
  
  int dcount = 0;
  int tcount = 0;
  int pcount = 0;
  int trcount = 0;
  int AdvancePhotonPointer;
  int DeleteMe, DeltaLevel, PauseMe;
  int prev_type = -1;
  float LightCrossingTime = RadiativeTransferRayMaximumLength * (VelocityUnits) /
    (clight * RadiativeTransferPropagationSpeedFraction); 
  FLOAT EndTime;
  if (RadiativeTransferAdaptiveTimestep)
    EndTime = PhotonTime+LightCrossingTime;
  else
    EndTime = PhotonTime+dtPhoton-PFLOAT_EPSILON;

  while (PP != NULL) {
    int retval = 0;
    if (PP->PreviousPackage == NULL)
      printf("Bad package.\n");
    DeleteMe = FALSE;
    PauseMe = FALSE;
    MoveToGrid = NULL;
    AdvancePhotonPointer = TRUE;
    if (MYPROC && DEBUG) {
      if(prev_type != PP->Type) {
	fprintf(stdout, "%s: Radiation type = %ld\n", __FUNCTION__, PP->Type);
	prev_type = PP->Type;
      }
    }
    if ((PP->CurrentTime) < EndTime) {
      retval = WalkPhotonPackage(&PP,
				 &MoveToGrid, ParentGrid, CurrentGrid, Grids0, nGrids0,
				 DeleteMe, PauseMe, DeltaLevel, LightCrossingTime,
				 LightSpeed, level, MinimumPhotonFlux);
      tcount++;
    } else {

      /* If all work is finished, store in FinishedPhotonPackages and
	 don't check for work until next timestep */

      SavedPP = PopPhoton(PP);
      PP = PP->NextPackage;
      InsertPhotonAfter(FPP, SavedPP);
      AdvancePhotonPointer = FALSE;

    }

    if (DEBUG > 1) 
      fprintf(stdout, "photon #%"ISYM" %x %x %x\n",
	      tcount,  PP,  PhotonPackages, 
	      MoveToGrid); 

    if (PauseMe == TRUE) {
      if (DEBUG > 1) fprintf(stdout, "paused photon %x\n", PP);
      this->RegridPausedPhotonPackage(&PP, ParentGrid, &MoveToGrid, DeltaLevel,
				      DeleteMe, DomainWidth, LightSpeed);

      // Insert in paused photon list if it belongs in this grid.
      if (MoveToGrid == NULL && DeleteMe == FALSE) {
	SavedPP = PopPhoton(PP);
	PP = PP->NextPackage;
	InsertPhotonAfter(PausedPP, SavedPP);
	AdvancePhotonPointer = FALSE;
      }
      pcount++;
    }

    if (DeleteMe == TRUE) {
      if (DEBUG > 1) fprintf(stdout, "delete photon %x\n", PP);
      dcount++;
      PP = DeletePhotonPackage(PP);
      MoveToGrid = NULL;
    } 

    if (MoveToGrid != NULL) {
      if (DEBUG > 1) {
	fprintf(stdout, "moving photon from %x to %x\n", 
		 CurrentGrid,  MoveToGrid);
	fprintf(stdout, "moving photon %x %x %x %x\n", 
		 PP,  PP->PreviousPackage, 
		 PP->NextPackage,  PhotonPackages);
      }
      ListOfPhotonsToMove *NewEntry = new ListOfPhotonsToMove;
      NewEntry->NextPackageToMove = (*PhotonsToMove)->NextPackageToMove;
      (*PhotonsToMove)->NextPackageToMove = NewEntry;
      NewEntry->PhotonPackage = PP;
      NewEntry->FromGrid = CurrentGrid;
      NewEntry->ToGrid   = MoveToGrid;
      NewEntry->ToGridNum= MoveToGrid->GetGridID();
      NewEntry->ToLevel  = level + DeltaLevel;
      NewEntry->ToProcessor = MoveToGrid->ReturnProcessorNumber();
      if (PauseMe)
	NewEntry->PausedPhoton = TRUE;
      else
	NewEntry->PausedPhoton = FALSE;
      if (NewEntry->ToProcessor >= NumberOfProcessors ||
	  NewEntry->ToProcessor < 0) {
	PP->PrintInfo();
	ENZO_VFAIL("Grid %d, Invalid ToProcessor P%d", GridNum, 
		   NewEntry->ToProcessor)
      }

      if (PP->PreviousPackage != NULL) 
	PP->PreviousPackage->NextPackage = PP->NextPackage;
      if (PP->NextPackage != NULL) 
	PP->NextPackage->PreviousPackage = PP->PreviousPackage;
      trcount++;
    } // ENDIF MoveToGrid

    if (AdvancePhotonPointer == TRUE)
      PP = PP->NextPackage;

  } // ENDWHILE photons

  if (DEBUG)
    fprintf(stdout, "grid::TransportPhotonPackage[%d]: "
	    "transported %"ISYM" deleted %"ISYM" paused %"ISYM" moved %"ISYM"\n",
	    this->ID, tcount, dcount, pcount, trcount);
  NumberOfPhotonPackages -= dcount;

#ifdef UNUSED
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    if (HasRadiation == TRUE) break;
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      if (HasRadiation == TRUE) break;
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	if (BaryonField[kphHINum][index] > 0) {

	  HasRadiation = TRUE;
	  break;
	}
      } // ENDFOR i
    }  // ENDFOR j
  } // ENDFOR k
#endif /* UNUSED */

  // Debug xyz-axis for a unigrid 64^3 with a source in the corner.
#define NO_DEBUG_AXES
#ifdef DEBUG_AXES
  printf("PHDebug(x): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][14914], BaryonField[kphHINum][14915], 
	 BaryonField[kphHINum][14916], 
	 BaryonField[kphHeIINum][14914], BaryonField[kphHeIINum][14915], 
	 BaryonField[kphHeIINum][14916], 
	 BaryonField[HINum][14914], BaryonField[HINum][14915], 
	 BaryonField[HINum][14916]);
  printf("PHDebug(y): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][14983], BaryonField[kphHINum][15053], 
	 BaryonField[kphHINum][15123], 
	 BaryonField[kphHeIINum][14983], BaryonField[kphHeIINum][15053], 
	 BaryonField[kphHeIINum][15123], 
	 BaryonField[HINum][14983], BaryonField[HINum][15053], 
	 BaryonField[HINum][15123]);
  printf("PHDebug(z): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][19813], BaryonField[kphHINum][24713], 
	 BaryonField[kphHINum][29613], 
	 BaryonField[kphHeIINum][19813], BaryonField[kphHeIINum][24713], 
	 BaryonField[kphHeIINum][29613], 
	 BaryonField[HINum][19813], BaryonField[HINum][24713], 
	 BaryonField[HINum][29613]);
#endif /* DEBUG_AXES */
	 
  return SUCCESS;
}
