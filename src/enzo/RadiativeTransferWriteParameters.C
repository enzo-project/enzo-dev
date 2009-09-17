/***********************************************************************
/
/  WRITES RADIATIVE TRANSFER PARAMETERS TO AN OUTPUT FILE
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* function prototypes */

int RadiativeTransferWriteParameters(FILE *fptr)
{

  fprintf(fptr, "RadiativeTransferRadiationPressure        = %"ISYM"\n", 
	  RadiationPressure);
  fprintf(fptr, "RadiativeTransferSourceRadius             = %"GSYM"\n", 
	  RadiativeTransferSourceRadius);
  fprintf(fptr, "RadiativeTransferPropagationSpeedFraction = %"GSYM"\n", 
	  RadiativeTransferPropagationSpeedFraction);
  fprintf(fptr, "RadiativeTransferPropagationDistance      = %"GOUTSYM"\n", 
	  RadiativeTransferPropagationDistance);
  fprintf(fptr, "RadiativeTransferCoupledRateSolver        = %"ISYM"\n", 
	  RadiativeTransferCoupledRateSolver);
  fprintf(fptr, "RadiativeTransferOpticallyThinH2          = %"ISYM"\n", 
	  RadiativeTransferOpticallyThinH2);
  fprintf(fptr, "RadiativeTransferPeriodicBoundary         = %"ISYM"\n", 
	  RadiativeTransferPeriodicBoundary);
  fprintf(fptr, "RadiativeTransferSplitPhotonRadius        = %"FSYM"\n", 
	  RadiativeTransferSplitPhotonRadius);
  fprintf(fptr, "RadiativeTransferRaysPerCell              = %"FSYM"\n", 
	  RadiativeTransferRaysPerCell);
  fprintf(fptr, "RadiativeTransferTimestepVelocityLimit    = %"FSYM"\n", 
	  RadiativeTransferTimestepVelocityLimit);
  fprintf(fptr, "RadiativeTransferInitialHEALPixLevel      = %"ISYM"\n", 
	  RadiativeTransferInitialHEALPixLevel);
  fprintf(fptr, "RadiativeTransferPhotonEscapeRadius       = %"FSYM"\n", 
	  RadiativeTransferPhotonEscapeRadius);
  fprintf(fptr, "RadiativeTransferInterpolateField         = %"ISYM"\n", 
	  RadiativeTransferInterpolateField);
  fprintf(fptr, "RadiativeTransferSourceClustering         = %"ISYM"\n", 
	  RadiativeTransferSourceClustering);
  fprintf(fptr, "RadiativeTransferPhotonMergeRadius        = %"FSYM"\n\n", 
	  RadiativeTransferPhotonMergeRadius);
  fprintf(fptr, "RadiativeTransferHIIRestrictedTimestep    = %"ISYM"\n", 
	  RadiativeTransferHIIRestrictedTimestep);
  
  return SUCCESS;
}
