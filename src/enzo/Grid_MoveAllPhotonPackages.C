#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (MOVE ALL PHOTONPACKAGES FROM SPECIFIED GRIDS TO THIS GRID)
/
/  written by: Tom Abel
/  date:       May, 2005
/  modified1:
/
/  PURPOSE: 
/
/    NOTE: We assume all the from grids are at the same level!
/    NOTE: Modelled after GBs Grid_MoveAllParticles.C
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);

int grid::MoveAllPhotonPackages(int NumberOfGrids, grid* FromGrid[])
{

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("grid::MoveAllPhotonPackages: NumberOfGrids(%"ISYM") must be > 0.\n", 
	    NumberOfGrids)
  }

  /* Determine total number of particles. */

  int TotalNumberOfPackages = NumberOfPhotonPackages;
  int i, j, gridcount, dim, *Number, *Type;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++) 
    TotalNumberOfPackages += FromGrid[gridcount]->NumberOfPhotonPackages;
  if (TotalNumberOfPackages == 0)
    return SUCCESS;

  /* Debugging info. */

//  if (debug)
//    fprintf(stdout, "MoveAllPackages: %"ISYM" (before: ThisGrid = %"ISYM").\n",
//	    TotalNumberOfPackages, NumberOfPhotonPackages);

  // go to end of List
  PhotonPackageEntry *PP = PhotonPackages->NextPackage;

  /* Error check number of photons.  If a bad value, reset photons */

  if (NumberOfPhotonPackages < 0) {
    printf("MoveAllPackages: WARNING. Resetting photons. "
	   "NumberOfPhotons = %"ISYM"\n", NumberOfPhotonPackages);
    NumberOfPhotonPackages = 0;
    TotalNumberOfPackages = 0;
    while (PP != NULL) {
      PP = DeletePhotonPackage(PP);
      PP = PP->NextPackage;
      TotalNumberOfPackages++;
    }
    printf("MoveAllPackages: deleted %"ISYM" photons\n", TotalNumberOfPackages);
    return SUCCESS;
  }

  /* Connect linked lists of FromGrids' PhotonPackages to this list */

  PP = PhotonPackages;

  int count = NumberOfPhotonPackages;
  int fromcount;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++) {

   /* If on the same processor, just add. */

    if (MyProcessorNumber == ProcessorNumber &&
        MyProcessorNumber == FromGrid[gridcount]->ProcessorNumber) {
      PhotonPackageEntry *FGPP = 
	(FromGrid[gridcount]->ReturnPhotonPackagePointer())->NextPackage;

      PhotonPackageEntry *NextFGPP = FGPP;

      fromcount = 0;

      while (FGPP != NULL) {

	// Get next package in from grid before transferring
	NextFGPP = FGPP->NextPackage;

	if (PP == NULL)
	  fprintf(stdout, "PP undefined.\n");
	FGPP->NextPackage = PP->NextPackage;
	if (PP->NextPackage != NULL)
	  PP->NextPackage->PreviousPackage = FGPP;
	PP->NextPackage = FGPP;
	FGPP->PreviousPackage = PP;

	FGPP = NextFGPP;

	fromcount++;
	count++;

	if (fromcount > FromGrid[gridcount]->ReturnNumberOfPhotonPackages()) {
	  printf("MoveAllPackages[P%"ISYM"]: WARNING! fromcount > #ph - %"ISYM" %"ISYM"\n",
		 MyProcessorNumber, fromcount,
		 FromGrid[gridcount]->ReturnNumberOfPhotonPackages());
	  printf("--> Ignoring the rest.\n");
	  NextFGPP = NULL;
	}

//	FGPP->PreviousPackage = PP;
//	PP->NextPackage = FGPP;
//	count++;
//	while (((PP->NextPackage)->NextPackage) != NULL) { 
//	  count++;
//	  PP=PP->NextPackage;
//	}

      } /* ENDWHILE FGPP != NULL */

      if (DEBUG)
	if (fromcount)
	  printf("MoveAllPackages[P%"ISYM"]: (LOCAL) counted %"ISYM" PhotonPackages. "
		 "grid #%"ISYM" of %"ISYM".\n", MyProcessorNumber, fromcount, gridcount, 
		 NumberOfGrids);

    }
    /* Otherwise, communicate. */
    
    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[gridcount]->ProcessorNumber) {
	if (DEBUG)
	  printf("MoveAllPackages(%"ISYM"): (COMM) moving %"ISYM" PhotonPackages. "
		 "grid #%"ISYM" of %"ISYM".\n", 
		 MyProcessorNumber, FromGrid[gridcount]->NumberOfPhotonPackages,
		 gridcount, NumberOfGrids);      
	if (FromGrid[gridcount]->CommunicationSendPhotonPackages(this, 
	       ProcessorNumber, NumberOfPhotonPackages, 
               FromGrid[gridcount]->NumberOfPhotonPackages, &PP) == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendPhotonPackages.\n");
	}
	count += FromGrid[gridcount]->ReturnNumberOfPhotonPackages();
	if (DEBUG)

	  printf("MoveAllPackages: (COMM) counted %"ISYM" PhotonPackages. "
		 "grid #%"ISYM" of %"ISYM".\n", count, gridcount, NumberOfGrids);      
      }
    }

  } // end: loop over grids.

  /* Set new number of particles in this grid. */

  NumberOfPhotonPackages = TotalNumberOfPackages;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++)
    FromGrid[gridcount]->SetNumberOfPhotonPackages(0);

  return SUCCESS;
}
