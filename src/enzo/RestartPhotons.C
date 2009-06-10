/***********************************************************************

  RESTART RADIATIVE TRANSFER
  written by:  John H. Wise
  date:        30 Mar 2006
  modified1:   
  
  Purpose: When we restart with radiation sources already in place, we
           must calculate the pre-existing radiation field before any
           other physics is cacluated.  We subtract a light-crossing
           time of the box from PhotonTime and run EvolvePhotons until
           the current time is reached.  This should populate the
           entire box with the correct radiation field.

 ***********************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
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
#include "StarParticleData.h"

#define CONVERGE 0.001

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars);

int RestartPhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		   Star *AllStars)
{

  int level;
  LevelHierarchyEntry *Temp;

  //MetaData->FirstTimestepAfterRestart = FALSE;

  if (AllStars == NULL)
    return SUCCESS;

  if (!RadiativeTransfer)
    return SUCCESS;

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    MassUnits, DensityUnits; 

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, MetaData->Time) == FAIL) {
    fprintf(stdout, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }
  
  /* Light crossing time */

  const float clight = 2.9979e10;

  float LightCrossingTime = VelocityUnits / 
    (clight * RadiativeTransferPropagationSpeedFraction);
  FLOAT SavedPhotonTime = PhotonTime;
  PhotonTime -= LightCrossingTime;

  if (debug)
    printf("Restarting radiative transfer.  Light-crossing time = %"GSYM"\n", 
	   LightCrossingTime);

  /* Solve radiative transfer */

  int PhotonCount, LastPhotonCount = 0;
  int savedCoupledChemistrySolver = RadiativeTransferCoupledRateSolver;
  RadiativeTransferCoupledRateSolver = FALSE;

  while ((dtPhoton > 0.) && RadiativeTransfer && 
	 (MetaData->Time >= PhotonTime))  {
    if (debug) 
      printf("EvolvePhotons[restart]: dt = %"GSYM", Time = %"FSYM", ", 
	     dtPhoton, PhotonTime);
    if (EvolvePhotons(MetaData, LevelArray, AllStars) == FAIL) {
      fprintf(stderr, "Error in EvolvePhotons.\n");
      ENZO_FAIL("");
    }

    PhotonCount = 0;
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->CountPhotonNumber();
	if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())
	  PhotonCount += Temp->GridData->ReturnNumberOfPhotonPackages();
	Temp = Temp->NextGridThisLevel;
      }
    }

#ifdef USE_MPI
    int value = PhotonCount;
    MPI_Allreduce(&value, &PhotonCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif /* USE_MPI */    

    if (LastPhotonCount > 0)
      if (float(PhotonCount-LastPhotonCount)/float(LastPhotonCount) < CONVERGE)
	PhotonTime = SavedPhotonTime + dtPhoton*1e-2;

    if (PhotonCount == 0 && LastPhotonCount == 0)
      PhotonTime = SavedPhotonTime + dtPhoton*1e-2;

    LastPhotonCount = PhotonCount;

  } /* ENDWHILE evolve photon */
  
  RadiativeTransferCoupledRateSolver = savedCoupledChemistrySolver;

  /* Optically thin Lyman-Werner (H2) radiation field */

  if (RadiativeTransferOpticallyThinH2)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->AddH2Dissociation(AllStars) == FAIL) {
	  fprintf(stderr, "Error in AddH2Dissociation.\n");
	  ENZO_FAIL("");
	}

  return SUCCESS;

}
